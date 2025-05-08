import pandas as pd


def read_input():
    #read input parameters
    input_dict={'num_muons':100,'random_seed':501}
    with open("input/input.txt") as f:
        for line in f:
            s=line.split()
            if (len(s)==2):
                if(s[0] in input_dict):
                    input_dict[s[0]]=int(s[1])

    #read energy loss file
    loss=pd.read_csv("input/loss.txt",sep=r'\s+',skiprows=26,header=None)
    loss.columns=['E', 'momentum', 'a', 'bb', 'bp','bph','br','rate','range','delta','beta','rate2']
    #change units
    loss['range']/=100.0
    loss['E']/=1000.0
    loss['rate']*=(100.0/1000)*1.03 # additional factor for the density of sea water
    return input_dict,loss