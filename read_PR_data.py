import numpy as np 
import matplotlib.pyplot as plt

def PR_2_dicts(dataFile,numRows,numCols):
    ''' 
    Each data (text) file can have several measurement types, so data for each measurement
    will be returned in its own dictionary. Also, depending on which plate reader was used, the data
    files are formatted slightly differently. So this code can handle both formats. 
    '''

    # generating the labels
    import string
    row_keys = list(string.ascii_uppercase[0:numRows])
    col_keys = range(1,numCols+1)

    from itertools import product
    cartprod = list(product(row_keys,col_keys))

    well_keys = []
    for i in cartprod:
        well_keys.append(str(i[0]+str(i[1])))

    dataRead = open(dataFile)

    count = 0
    ReadsLine = [] # get the line in which the number of reads appears. 
    TimeList = [] # where 'Time' appears in the data file is where data starts (except the first appearance), and where 'Results' appears is where the data ends. 
    for line in dataRead:
        if 'Reads' in line: 
            ReadsLine.append(line)
        if 'Time' in line:
            TimeList.append(count) 
        if 'Results' in line:
            TimeList.append(count) 
            break
        count += 1

    Reads = [] # getting the number of reads
    for word in ReadsLine[0].split():
        if word.isdigit():
            Reads.append(int(word))

    dataFinal = []  
    # need to properly define where the different measurements start and end
    for k in range(1,len(TimeList)-1):

        if k == (len(TimeList) - 2):
            dataStart = TimeList[k]
            dataEnd = TimeList[k+1] - 1
        else:
            dataStart = TimeList[k]
            dataEnd = TimeList[k+1] - 4

        dataRead = open(dataFile) 

        # creating a list of the data
        count = 0 
        data = []
        for line in dataRead:
            if count >= dataStart and count <= dataEnd:
                data.append(line)
            count += 1

        # formatting the data appropriately
        for i in range(len(data)):
            data[i] = data[i].split('\t') # new element in each line after tab
            data[i][-1] = data[i][-1].strip() # get rid of '\n' at end of each line
            data[i] = data[i][2:] # don't want the time and temp elements
            data[i] = ' '.join(data[i]).split() # remove the empty strings 

        # removing all empty lists
        data = [x for x in data if x != []]

        # replacing any 'OVRFLW' elements with 0.0
        for ii in range(len(data)):
            for jj in range(len(data[ii])): 
                if data[ii][jj] == 'OVRFLW':
                    data[ii][jj] = 0.0

        # convert to floats and don't need the first list that contains all the well labels
        data = np.array(data[1:],dtype=np.float) # each row contains a snapshot for some time, and each column contains the trajectory of a single state

        # assign the key-value pairs to a dictionary
        data_dict = {}
        for i in range(len(well_keys)):
            data_dict[well_keys[i]] = data[:,i:i+1]


        dataFinal.append(data_dict)

    return dataFinal,Reads[0]

# # Uncomment below for an example. 


# # The only manual inputs that are necessary are the number of columns, number of rows
# # For plotting purposes, also supply the sampling rate. 

# dataFile = 'pf_wt_growth_8_7_AH.txt'
# numRows = 2
# numCols = 5
# sampling_rate = 5 # minutes
# data,numReads = PR_2_dicts(dataFile,numRows,numCols) 
# ''' data will contain a list of arrays, length of 
#     this list will be length of how many data types you measured 
#     in the plate reader e.g. if you measured OD and fluorescence, then data will be of length two.
#     Whichever measurement type comes first in the data file will come first in the list. '''
# data_od600 = data[0]
# # data_fluor = data[1] # if there was another measurement taken on the plate reader (e.g. fluorescence)
# tSpan = np.linspace(0,numReads*sampling_rate,numReads)/60 # represent time in hours


# # # plot some wells in row A of plate by specifying dictionary keys:

# plt.figure();
# plt.plot(tSpan,data_od600['A1']);
# plt.plot(tSpan,data_od600['A3']);
# plt.plot(tSpan,data_od600['A5']);
# plt.show();


# # plot all wells by generating the dictionary keys and placing them in a list: 

# import string
# row_keys = list(string.ascii_uppercase[0:numRows])
# col_keys = range(1,numCols+1)

# from itertools import product
# cartprod = list(product(row_keys,col_keys))

# well_keys = []
# for i in cartprod:
#   well_keys.append(str(i[0]+str(i[1])))

# # and then you can loop through the keys or a subset of them for more concise plotting

# plt.figure();
# for well in well_keys: 
#   plt.plot(tSpan,data_od600[well])
# plt.show();





