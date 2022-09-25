from time import perf_counter
from numpy import log10, array, where

def convert_time(time_second):
    exponent = round(log10(time_second))
    exp_arr = -array([0, 3, 6, 9, 12])
    exp_name = {0:'', -3:'m', -6:'µ', -9:'n', -12:'p'}
    exp_close = abs(exp_arr - exponent)
    exp_index = where(exp_close == min(exp_close))
    exp_fix = exp_arr[exp_index][0]
    exp_fix10 = 10. ** exp_fix
    return f'{time_second / exp_fix10:8.5f} {exp_name[exp_fix]}s'

def timer(func):
    def inner(*args, **kwargs):
        start = perf_counter()
        func_result = func(*args, **kwargs)
        end = perf_counter()
        time_diff = end - start
        print(f'{func.__name__}: {convert_time(time_diff)}')
        return func_result
    return inner


def time_it(run_num = 10, loop_num = 1000):
    def inner(func):
        def wrapper(*args, **kwargs):
            avg_times = []
            for i in range(run_num):
                start = perf_counter()
                for j in range(loop_num):
                    func_result = func(*args, **kwargs)
                end = perf_counter()
                avg_times += [(end - start)/loop_num]
            avg_time = sum(avg_times)/run_num
            std_dev = sum([(t-avg_time)**2 for t in avg_times])/(run_num-1)
            cavg_time = convert_time(avg_time)
            cstd_dev = convert_time(std_dev)
            print(f'{func.__name__}: run_num = {run_num}, loop_num = {loop_num}, avg time = {cavg_time} ± {cstd_dev}')
            return func_result
        return wrapper
    return inner

def timer_str(f):
    start = perf_counter()
    result = eval(f)
    end = perf_counter()
    print(end - start)
    return result

if __name__=='__main__':
    @time_it()
    def for_loop():
        x = 0
        for i in range(10):
            x += 1
            x *= 2
            x /= 3
        return x

    @timer
    def for_loop2():
        x = 0
        for i in range(10):
            x += 1
            x *= 2
            x /= 3
        return x

    for i in range(10):
        print(i, end=' ')
        for_loop()

    for i in range(10):
        print(i, end=' ')
        for_loop2()