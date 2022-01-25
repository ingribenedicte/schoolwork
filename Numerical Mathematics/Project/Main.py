import Functions as func
import Submodules

def main():
    method = input("Name one of the following methods: 'Newton', 'Improved Newton', 'Olver', 'Improved Olver'")
    function = input("Name one of the functions: 'function 1', function 2', function 3'")
    if function == 'function 1':
        f = func.f1
        df = func.df1
        ddf = func.ddf1
    if function == 'function 1':
        f = func.f2
        df = func.df2
        ddf = func.ddf2
    if function == 'function 2':
        f = func.f2
        df = func.df2
        ddf = func.ddf2
    if function == 'function 3':
        f = f3
        df = df3
        ddf = ddf3