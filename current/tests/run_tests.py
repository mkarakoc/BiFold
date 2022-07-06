import os

def print_title(title='', n_l='\n', omit=''):
    l_title = len(title)
    sides = (110-l_title)//2
    if sides>1:
        side_strip = '═' * (sides-1)
        print(f'{omit}{n_l}{omit}{side_strip}╣{title}╠{side_strip}{omit}{n_l}{omit}')
    else:
        print(f'{omit}{n_l}{omit}{title}')

working_dir = os.getcwd()
for dirpath, dnames, fnames in os.walk("./"):
    for fname in fnames:
        if fname.endswith(".py") and 'run_test' not in fname:
            test_file = os.path.join(dirpath, fname)
            print_title(title=test_file)
            os.chdir(dirpath)
            os.system(f'python {fname}')
            os.chdir(working_dir)
