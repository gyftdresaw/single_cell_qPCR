
# results data in '20140217 PLZF Fetal.xlsx' (not normalized to housekeeping genes)
# is currently split into the results from one 96.96 chip and the combined results
# of two 48.48 chips

# This script will:
# => un-combine the results from the two smaller chips inorder to maybe 
#    see how comparable results are from different chips
# => combine all the chip results into one big matrix for convenience

import numpy as np
import string

## splitting results from chips 2 and 3 ##
chip23_datafile = 'data/parsed/chip2-3_results.txt'

# file lines are separated by carriage returns '\r'
# need to open file with universal newlines enabled (everything treated as '\n')
infile = open(chip23_datafile,'rU')

data23 = [line.strip().split('\t') for line in infile.readlines()]
data23_array = np.array(data23) # for easier indexing
data23_array[data23_array == 'N.D.'] = 'NA' # replacing with NAs
r,c = data23_array.shape

# get indices associated with each chip
chip2_indices = [i+j*12 for j in xrange(8) for i in xrange(1,7)]
chip3_indices = [i+6+j*12 for j in xrange(8) for i in xrange(1,7)]

# split up data according to chip
data2_array = np.copy(data23_array[:,[0]+chip2_indices])
data3_array = np.copy(data23_array[:,[0]+chip3_indices])

# write these to file
# format is '%s' and delimiter is '\t'
np.savetxt('data/parsed/chip2_clean.txt',data2_array[[0]+range(2,r),:],'%s','\t')
np.savetxt('data/parsed/chip3_clean.txt',data3_array[[0]+range(2,r),:],'%s','\t')

infile.close()

## obtain chip1 data ##
chip1_datafile = 'data/parsed/chip1_results.txt'

infile = open(chip1_datafile,'rU')
data1 = [line.strip().split('\t') for line in infile.readlines()]
data1_array = np.array(data1) # easier indexing
data1_array[data1_array == 'N.D.'] = 'NA' # replacing with NAs
r1,c1 = data1_array.shape

# write this to file
np.savetxt('data/parsed/chip1_clean.txt',data1_array,'%s','\t')

infile.close()

## combine all data into one big array ##
combined_array = np.concatenate((data1_array,data23_array[3:,:]))
np.savetxt('data/parsed/combined_clean.txt',combined_array,'%s','\t')
