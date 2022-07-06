import os
from bifold import print_title

workdir = os.getcwd()
for dirpath, dnames, fnames in os.walk("./"):
    for fname in fnames:
        if fname.endswith(".py") and 'run_test' not in fname:
            test_file = os.path.join(dirpath, fname)
            print_title(title=test_file)
            os.chdir(dirpath)
            os.system(f'python {fname}')
            os.chdir(workdir)
