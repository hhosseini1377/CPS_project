from fpEDF_VD import generate_tasks
import matplotlib.pylab as plt

def fpEDF_VD_rho():
    """
        Compare results of the fpEDF_VD algorithm in response to the rho 
    """
    work_loads = generate_tasks()
    total_tasks_per_rho = {
        0.3: 0,
        0.5: 0,
        0.7: 0,
        0.9: 0,
    }

    total_schedulable_tasks_per_rho = {
        0.3: 0,
        0.5: 0,
        0.7: 0,
        0.9: 0,
    }
    
    for task_set in work_loads:
        task_set.fpEDF_VD_preprocessing()
        total_tasks_per_rho[task_set.rho] += 1
        if task_set.schedulable:
            total_schedulable_tasks_per_rho[task_set.rho] += 1
    
    schedulibility_ratio_rho = {}
    for rho in total_tasks_per_rho:
        schedulibility_ratio_rho[rho] = total_schedulable_tasks_per_rho[rho]/total_tasks_per_rho[rho]
    plt.bar(range(len(schedulibility_ratio_rho)), list(schedulibility_ratio_rho.values()), align='center')
    plt.xticks(range(len(schedulibility_ratio_rho)), list(schedulibility_ratio_rho.keys()))
    plt.show()

def fpEDF_VD_utilization():
    work_loads = generate_tasks()
    total_tasks_per_utilization = {}
    total_schedulable_tasks_per_utilization = {}

    range = 0 

    while range < 1:
        total_tasks_per_utilization[range] = 0
        total_schedulable_tasks_per_utilization[range] = 0
        range += 0.05
        range = round(range, 2)

    for task_set in work_loads:
        task_set.fpEDF_VD_preprocessing()
        total_tasks_per_utilization[round((task_set.utilization[0]//0.05)* 0.05, 2)] += 1
        if task_set.schedulable:  
            total_schedulable_tasks_per_utilization[round((task_set.utilization[0]//0.05)* 0.05, 2)] += 1

    schedulability_ratio_utilization = {}
    for utilization in  total_tasks_per_utilization:
        if total_tasks_per_utilization[utilization] != 0:
            schedulability_ratio_utilization[utilization] = total_schedulable_tasks_per_utilization[utilization]/total_tasks_per_utilization[utilization]
        else:
            schedulability_ratio_utilization[utilization] = None

    print(schedulability_ratio_utilization)
    print(list(schedulability_ratio_utilization.values()))
    util, ratio = zip(*schedulability_ratio_utilization.items())
    plt.bar(util, ratio)
    plt.show()
    # plt.bar(range(len(schedulability_ratio_utilization)), list(schedulability_ratio_utilization.values()), align='center')
    # plt.xticks(range(len(schedulability_ratio_utilization)), list(schedulability_ratio_utilization.keys()))
    # plt.show()

fpEDF_VD_utilization()