from UUnifast import UUniFastDiscard
from numpy import random 
import random as rand
import numpy as np
import cvxpy as cvx

# Define constant values
n_set = [20, 40, 60, 80, 100]
m_set = [4]
rho_set = [0.6]
C_down, C_up = 1, 100
R, P = 4, 0.5
task_sets_num = 1000
time_delta = 0.01
time_decimals = 2

# Define a class for tasks
class task():

    def __init__(self, u_h) -> None:

        # Determine the criticality of the task
        criticality_rand = random.uniform(0, 1, 1)[0]
        if criticality_rand < P:
            self.high_critical = True
        else:
            self.high_critical = False
            
        self.u_h = u_h

        # Determine the attributes if the criticality of the task is low 
        if not self.high_critical:
            self.u_l = u_h
            self.C_l = round(random.uniform(C_down, C_up, 1)[0], time_decimals)
            self.C_h = round(self.C_l, time_decimals) 
            self.period = round(self.C_h/u_h, time_decimals)

        # Determine the attributes if the criticality of the task is high
        else:
            self.u_l = random.uniform(u_h, u_h/R, 1)[0]
            self.C_l = round(random.uniform(C_down, C_up, 1)[0], time_decimals)
            self.period = round(self.C_l/self.u_l, time_decimals)
            self.C_h = round(self.period*self.u_h, time_decimals)
        self.next_release = self.period
    
    def set_virtual_deadline(self, x):
        if self.high_critical:
            self.virtual_deadline = round(self.period*x, time_decimals)
        else:
            self.virtual_deadline = self.period
    
    def check_released(self, time):
        if time >= self.next_release:
            self.next_release += self.period
            return True
        return False


     
class task_job():
    def __init__(self, task, arrival_time, x) -> None:
        self.task = task
        self.arrival_time = arrival_time
        self.deadline = arrival_time + task.period
        self.virtual_deadline = arrival_time + round(task.period*x, time_decimals)
        self.completed = False
        self.deadline_missed = False
        self.preemptable = True

        # The actual execution time of a job
        # self.remaining_execution_time = round(random.uniform(0, task.C_h, 1)[0], time_decimals)
        self.remaining_execution_time = round(random.uniform(task.C_l/2, task.C_l, 1)[0], time_decimals)
        self.execution_time = 0

    """
    This function reduces remaining execution time of the job and 
    checks if the task is completed.
    """
    def check_completion(self, rho, time_delta, system_mode):
        if system_mode == 'high_criticality':
            speed_rate = 1
        else:
            speed_rate = rho
        self.remaining_execution_time -= round(time_delta*speed_rate, time_decimals)
        self.remaining_execution_time = round(self.remaining_execution_time, time_decimals)
        self.execution_time += round(time_delta*speed_rate, time_decimals)
        self.execution_time = round(self.execution_time, time_decimals)
        if self.remaining_execution_time <= 0:
            self.completed = True
            return True
        return False
    
    def check_low_criticality_miss(self):
        if self.execution_time > self.task.C_l:
            return True
        else:
            return False
    
    def check_high_criticality_miss(self, time):
        if time > self.deadline:
            return True
        else:
            return False
    
    def check_running_on_core(self, cores, ready_queue, system_mode):
        # Check if any core is empty
        for core_ in cores:
            if not core_.running_job:
                core_.running_job = self
                ready_queue.remove(self)
                return True

        replace_with_highest_deadline_job(ready_queue=ready_queue, system_mode=system_mode, job=self, cores=cores, is_laxity=False)
        return False
            
# Define a class for task sets
class work_load():

    def __init__(self, n, m, utilization, rho) -> None:
        self.set = []
        self.n = n
        self.m = m
        self.utilization = utilization
        self.rho = rho
    
    """
    This function determines whether the workload is schedulable or not.
    if it is schedulable, the ratio for the virtual deadline will be set.
    """
    def fpEDF_VD_preprocessing(self):

        tasks = self.set
        m = self.m
        uLOvalues = np.array([t.u_l for t in tasks])
        uHIvalues = np.array([t.u_h for t in tasks])
        uL_total = sum(uLOvalues)
        uH_total = sum(uHIvalues)
        uL_max = max(uLOvalues)
        uH_max = max(uHIvalues)

        termA = max(uL_max/self.rho , uL_total/((self.m+1)*self.rho/2))
        termB = max(uH_max, uH_total/((self.m+1)/2))

        self.x = termA
        output_fpEDFVD = ((termA + termB) <= 1)

        return output_fpEDFVD

    def MCF_FR_pre_processing(self):

        tasks = self.set
        m = self.m
        uLOvalues = np.array([t.u_l for t in tasks])
        uHIvalues = np.array([t.u_h for t in tasks])
        uL_total = sum(uLOvalues)
        uH_total = sum(uHIvalues)
        uL_max = max(uLOvalues)
        uH_max = max(uHIvalues)

        termA = uL_total/(m + uL_total - uH_total)
        termB = max(uLOvalues/(1+uLOvalues-uHIvalues))
        lambdaValue = max(termA, termB)

        if uL_total >= self.rho*self.m:
            output_MCFFR = False
        elif uH_total >= m:
            output_MCFFR = False
        elif lambdaValue <= self.rho:
            output_MCFFR = True
            theta3 = uLOvalues/lambdaValue + uHIvalues - uLOvalues
            theta_L3 = lambdaValue * theta3
            theta_H3 = theta3
        else:
            output_MCFFR = False

        return output_MCFFR
    
    def MCF_MP_pre_processing(self, ):
        n = self.n
        tasks = self.set

        uLOvalues = np.array([t.u_l for t in tasks])
        uHIvalues = np.array([t.u_h for t in tasks])
        cLOvalues =  np.array([t.C_l for t in tasks])
        cHIvalues =  np.array([t.C_h for t in tasks])
        tvalues =  np.array([t.period for t in tasks])

        # Creating CVX variables
        variables = list()
        # Creating CVX inverse variables for the non-linear (convex) condition
        invVariables = list()

        # Creating theta_LO variables
        for i in range(n):
            a = cvx.Variable()
            variables.append(a)
            invVariables.append(cvx.inv_pos(a))
        
        # Creating theta_HI variables
        for i in range(n):
            b = cvx.Variable()
            variables.append(b)
            invVariables.append(cvx.inv_pos(b))

        # Creating constraints
        constraints = list()
        for i in range(n):
            constraints.append(variables[i] <= self.rho)
            constraints.append(variables[i] >= uLOvalues[i])
            constraints.append(variables[i+n] <= 1)
            constraints.append(variables[i+n] >= uHIvalues[i])
            constraints.append(variables[i] <= variables[i+n])
            constraints.append(cLOvalues[i]*invVariables[i] + (cHIvalues[i] - cLOvalues[i])*invVariables[i+n] <= tvalues[i])
        constraints.append(cvx.sum([variables[i] for i in range(n)]) <= self.rho*self.m)
        constraints.append(cvx.sum([variables[i+n] for i in range(n)]) <= self.m)

        # Creating dummy objective function
        objective = cvx.Minimize(1)

        # CVX problem formulation
        CVX_FAILURE = False
        problem = cvx.Problem(objective, constraints)
        if problem.is_dcp():
            try:
                problem.solve(solver=cvx.CVXOPT)
            except cvx.error.SolverError:
                try:
                    problem.solve()
                except:
                    CVX_FAILURE = True
        else:
            CVX_FAILURE = True

        # if np.isinf(problem.value):
        if problem.status in ["infeasible", "unbounded"]:
            # No feasible set
            output_MCFMP_CVXPY = False
        else:
            output_MCFMP_CVXPY = True
        if not CVX_FAILURE:
            return output_MCFMP_CVXPY
        else:
            return 'error'


class core():
    def __init__(self) -> None:
        self.running_job = None

    # Remove the current job from the active jobs and change the current job
    def change_job_low_critical(self, ready_queue, active_queue):
        job = find_least_virtual_deadline(ready_queue)
        if self.running_job:
            active_queue.remove(self.running_job)
        self.running_job = job

    def change_job_high_critical(self, ready_queue, active_queue):
        job = find_least_actual_deadline(ready_queue)
        if self.running_job:
            active_queue.remove(self.running_job)
        self.running_job = job
            

# Generate the tasks
def generate_tasks():

    sets = []
    while len(sets) < task_sets_num:

        m = rand.choice(m_set)
        n = rand.choice(n_set)
        rho = rand.choice(rho_set)
        per_core_utilization = random.uniform(0, 1, 1)

        set = UUniFastDiscard(n, per_core_utilization, 1)[0]
        utilization_exceeded = False

        # Check if each utilization is less than 1. If this condition is met the set will be added to the task_sets.
        for i in range(len(set)):
            set[i] *= m
            if set[i] > 1:
                # Discard a task if its utilization exceeds 1.
                utilization_exceeded = True
                continue

        if not utilization_exceeded:
            generated_take_set = work_load(n, m, per_core_utilization, rho, )
            for util in set:
                generated_take_set.set.append(task(util[0]))
            sets.append(generated_take_set)
    print(sum(set))
    return sets

def find_least_virtual_deadline(ready_queue):
    """
    Return the job in the ready queue with the least virtual deadline. 
    """
    if ready_queue:
        least_deadline_job = ready_queue[0]
        for job in ready_queue:
            if job.virtual_deadline < least_deadline_job.virtual_deadline:
                least_deadline_job = job
        ready_queue.remove(least_deadline_job)
        return least_deadline_job
    return None

def replace_with_highest_deadline_job(ready_queue, system_mode, cores, job, is_laxity):
    """
    Find the preemptive task with highest deadline to be replaced by a job.
    """ 
    most_deadline = 0
    most_deadline_running_job = None  
    if system_mode == 'low_criticality':
        for core_ in cores:
            if core_.running_job and core_.running_job.preemptable:
                if core_.running_job.virtual_deadline > most_deadline:
                    most_deadline_running_job = core_.running_job
                    most_deadline = core_.running_job.virtual_deadline
        if most_deadline_running_job:
            if (job.virtual_deadline < most_deadline and not is_laxity) or is_laxity:
                ready_queue.append(core_.running_job)
                ready_queue.remove(job)
                core_.running_job = job
                return True
    else:
        for core_ in cores:
            if core_.running_job and core_.running_job.preemptable:
                if core_.running_job.deadline > most_deadline:
                    most_deadline_running_job = core_.running_job
                    most_deadline = core_.running_job.deadline
        if most_deadline_running_job:
            if (job.deadline < most_deadline and not is_laxity) or is_laxity:
                ready_queue.append(core_.running_job)
                ready_queue.remove(job)
                core_.running_job = job
                return True  
    
def find_least_actual_deadline(ready_queue):
    """
    Return the job in the ready queue with the least actual deadline.
    """
    if ready_queue:
        least_deadline_job = ready_queue[0]
        for job in ready_queue:
            if job.deadline < least_deadline_job.deadline:
                least_deadline_job = job
        ready_queue.remove(least_deadline_job)
        return least_deadline_job
    return None

def change_mode(active_queue, cores_list, destination_mode):
    ready_queue = active_queue.copy()

    # Remove jobs from  cores
    for core_ in cores_list:
        core_.running_job = None
    # Make the jobs preemptable
    for job in active_queue:
        job.preemptable = True

    # Remap jobs on  cores
    for core_ in cores_list:
        if destination_mode == 'low_criticality':
            least_deadline = find_least_virtual_deadline(ready_queue=ready_queue)
        else:
            least_deadline = find_least_actual_deadline(ready_queue=ready_queue)
        if least_deadline:
            core_.running_job = least_deadline
        else:
            break    
    return ready_queue

def check_laxity(ready_queue, time, rho, time_delta, system_mode, cores):
    
    for job in ready_queue:
        if system_mode == 'low_criticality':
            if job.remaining_execution_time >= time_delta + (job.virtual_deadline - time)/rho:
                replace_with_highest_deadline_job(ready_queue=ready_queue, system_mode=system_mode, cores=cores, job=job, is_laxity=True)
                job.preemptable = False
                print('There is no enough laxity, job set on the core')
        else:
            if job.remaining_execution_time >= time_delta + (job.deadline - time):
                replace_with_highest_deadline_job(ready_queue=ready_queue, system_mode=system_mode, cores=cores, job=job, is_laxity=True)
                job.preemptable = False 
                print('There is no enough laxity, job set on the core')

def fpEDF_VD(work_load):

    cores_mode = 'low_criticality'
    cores_num, rho = work_load.m, work_load.rho
    ready_queue = []
    active_queue = []
    tasks_queue = []
    cores_list = []

    # Initialize cores
    for _ in range(cores_num):
        new_core = core()
        cores_list.append(new_core)

    # Initialize jobs
    for task in work_load.set:
        new_job = task_job(task=task, arrival_time=0, x=work_load.x)
        ready_queue.append(new_job)
        active_queue.append(new_job)
        tasks_queue.append(new_job)
    
    # Initialize jobs on cores
    for core_ in cores_list:
        if len(ready_queue) > 0:
            # Select the job with least virtual deadline
            selected_job = find_least_virtual_deadline(ready_queue=ready_queue)
            core_.running_job = selected_job
    
    time = 0
    work_load_missed = False
    time_interval = 0

    while True:
        len(active_queue)
        check_laxity(ready_queue=ready_queue, time=time, rho=rho, time_delta=time_delta, system_mode=cores_mode, cores=cores_list)

        # Check if any task releases new job
        for task in work_load.set:
            if task.check_released(time):
                job = task_job(task=task, arrival_time=time, x=work_load.x)
                ready_queue.append(job)
                active_queue.append(job)
                tasks_queue.append(job)
                job.check_running_on_core(cores= cores_list, ready_queue=ready_queue, system_mode=cores_mode)


        if cores_mode == 'low_criticality':

            # Check if any task has completed its execution time and change the job running on the core (if exists)
            for core_ in cores_list:
                if core_.running_job:
                    if core_.running_job.check_completion(rho=rho, time_delta=time_delta, system_mode=cores_mode):
                        core_.change_job_low_critical(ready_queue, active_queue, )

            # Check if any task has missed its virtual deadline
            for job in active_queue:
                if job.check_low_criticality_miss():
                    if not job.task.high_critical:
                        work_load_missed = True
                    else:
                        print("System mode swtiched to HIGH CRITICALITY")
                        cores_mode = 'high_criticality'
                        ready_queue = change_mode(destination_mode=cores_mode, active_queue=active_queue, cores_list=cores_list)
                    
        elif cores_mode == 'high_criticality':

            # Check if any task has completed its execution time and change the job running on the core (if exists)
            for core_ in cores_list:
                if core_.running_job:
                    if core_.running_job.check_completion(rho=rho, time_delta=time_delta, system_mode=cores_mode):
                        core_.change_job_high_critical(ready_queue, active_queue, )

            # Check if any task has missed its actual deadline
            for job in active_queue:
                if job.check_high_criticality_miss(time):
                    work_load_missed = True
            
            # Check if no task exists in the system. If this condition is met change systems mode to low ciriticality
            if not active_queue:
                cores_mode = 'low_criticality'
                ready_queue = change_mode(destination_mode=cores_mode, active_queue=active_queue, cores_list=cores_list)
                print('System mode swtiched to LOW CRITICALITY')


        if work_load_missed:
            print('System missed its deadline at time', time)
            break
            
        time += round(time_delta, time_decimals)
        time = round(time, time_decimals)

        if time == time_interval:
            print(time)
            time_interval += 1


work_loads = generate_tasks()
results = []

i = 0
while i < 1:
    i = round(i, 2)
    result = {
        'utilization': i,
        'fpEDF-VD': 0,
        'MCF-FR': 0,
        'MCF-MP': 0,
        'count': 0
    }
    i += 0.05
    results.append(result)

total_task_set = 0
for task_set in work_loads:
    print(total_task_set)
    task_set_utilization = (task_set.utilization[0]//0.05)*0.05
    task_set_utilization = round(task_set_utilization, 2)
    
    for result in results:
        if result['utilization'] == task_set_utilization:
            break

    out_MCF_mp = task_set.MCF_MP_pre_processing()
    if out_MCF_mp != 'error':
        if out_MCF_mp:
            result['MCF-MP'] += 1
    else:
        continue
    if task_set.MCF_FR_pre_processing():
        result['MCF-FR'] += 1
    if task_set.fpEDF_VD_preprocessing():
        result['fpEDF-VD'] += 1

    total_task_set += 1
    result['count'] += 1

for result in results:
    result['MCF-FR'] = result['MCF-FR']  / result['count']
    result['MCF-MP'] = result['MCF-MP']  / result['count']
    result['fpEDF-VD'] = result['fpEDF-VD']  / result['count']
print(results)