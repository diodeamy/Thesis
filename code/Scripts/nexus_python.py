#!/usr/bin/python


#################################################################################################################################
#
# Execute the commands necessary to compute the NEXUS+ environments. The steps are the following:
#           0 - compute the maximum node signature
#           1 - compute the optimal detection threshold for nodes and write the valid node regions to a 'clean' file
#           2 - compute the maximum filament and wall signature
#           3 - compute the optimal detection threshold for filaments and write the valid filament regions to a 'clean' file
#           4 - compute the optimal detection threshold for walls and write the valid wall regions to a 'clean' file
#           5 - combine the node, filament and wall 'clean' files to a combined file with: 0=voids, 2=walls, 3=filaments and 4=nodes
# 
# Calling the script via 'nexus_script.py' will execute all the above steps in the given order.
# You can call the script with 2 additional options: 'nexus_script.py N1 N2' to execute steps N1 to N2 (including).
# Therefore to run steps 1 to 4 use: 'nexus_script.py 1 4'
#
#################################################################################################################################



# the path to the directory with the executable files
binDir    = 'Users/users/nastase/Applications/NEXUS_python/NEXUS_1.0/bin/'
# the input density file
densityFile = 'density.a_den'
# uncomment the following if to write the console messages to a log file. You can select a different log file name.
#logFile = 'output.log'



# Once you have the density, you can run the NEXUS algorithm
# the following computes the node/cluster maximum signature using the NEXUS_den method. We don't use NEXUS++ since it does not fare that well for cluster detection. The option '--den' selects the NEXUS_den method.
# The program takes as input the density file and a range of scales that give the set of scales used in the multiscale analysis. In this case, the set of scales '0-6' tell the program to compute the node signature for smoothing scales R=R_0 sqqrt(2)^N with N=0 to 6. These 7 smoothing scales are used to compute the overall scale independent signature. Using the default arguments in the program, the filter scales range from 0.5 to 4 Mpc/h.
# The option '--node' tells the code to compute the signature for nodes while the next argument gives the root name for the output file. In this case the output file is 'node_maxResponse.MMF'.
nexus_response_1 = '{0}/MMF_response {1} 0-6 --den --node node'.format( binDir, densityFile )


# Once we have the node maximum signature, we run 'MMF_threshold' to get the optimal threshold for node detection
# 'node_maxResponse.MMF' is the input maximum signature file while 'node.threshold' is an output text file containing the data used to get the optimal threshold. The value of the optimal threshold using the method described in the NEXUS paper is given in the 1st line of file  'node.threshold'.
# Option '--cleanFile node_clean.MMF' tells the program to output a file with all the regions above the threshold as values of 1 - this corresponds to valid clusters while the remaining voxels have a value of 0.
# Option '--node' specify that the threshold should be computed for nodes.
# Option '--virDen 350' gives the value of the virial density with respect to background density at the snapshot redshift. 
# Option '--minSize 5.e13' gives the smallest mass (in Msolar/h units) that an object must have to be considered a valid cluster.
# Option '--densityFile {1}' gives the name of the density file used to obtain the node signature. This is used to compute the mass inside cluster regions.
nexus_threshold_1 = '{0}/MMF_threshold node_maxResponse.MMF node.threshold --node --minSize 5.e13 --virDen 350 --cleanFile node_clean.MMF --densityFile {1}'.format( binDir, densityFile )


# This step computes the filament and wall maximum response using the NEXUS+ method (program option '--logFilter').
# The option '--fila fila' specifies to compute the filament maximum signature and to output it to the file 'fila_maxResponse.MMF'.
# The option '--wall wall' specifies to compute the wall maximum signature and to output it to the file 'wall_maxResponse.MMF'.
# For more details see the comments for command: nexus_response_1
nexus_response_2 = '{0}/MMF_response {1} 0-6 --logFilter --fila fila --wall wall'.format( binDir, densityFile )


# This steps computes the optimal detection threshold for filaments
# Takes as input the filament maximum signature 'fila_maxResponse.MMF' and outputs the data used for threshold computation in file 'fila.threshold'.
# Option '--fila' specifies that the input file is a filament signature one.
# Option '--minSize 10' specify the minimum volume in (Mpc/h)^3 that a distinct filament must have to be considered a valid feature. Distinct objects with smaller volume are discarded.
# Option '--cleanFile fila_clean.MMF' specifies to output the clean filament file. The file contains a value of 1 for filament voxels and 0 otherwise.
# Option '--densityFile {1}' gives the density file and is used to compute the mass in filaments and how it changes with filament signature.
# Option '--nodeFile node_clean.MMF' gives the node clean file compute in a previous step that is used to mask the node regions. This means that a cluster region will not be identified also as a filament.
nexus_threshold_2 = '{0}/MMF_threshold fila_maxResponse.MMF fila.threshold --fila --minSize 10 --cleanFile fila_clean.MMF --densityFile {1} --nodeFile node_clean.MMF'.format( binDir, densityFile )


# This steps computes the optimal detection threshold for walls.
# Takes as input the wall maximum signature 'wall_maxResponse.MMF' and outputs the data used for threshold computation in file 'wall.threshold'.
# The options are exactly the same as for the filament detection thershold, so see description above.
# Option '--filaFile fila_clean.MMF' gives the fila clean file compute in a previous step that is used to mask the filamentary regions. This means that a filament region will not be identified also as a wall region.
nexus_threshold_3 = '{0}/MMF_threshold wall_maxResponse.MMF wall.threshold --wall --minSize 10 --cleanFile wall_clean.MMF --densityFile {1} --nodeFile node_clean.MMF --filaFile fila_clean.MMF '.format( binDir, densityFile )


# This last steps combines the 3 clean files for nodes, filaments and walls in a single environment file.
# Takes as input options: node clean file, filament clean file and wall clean file.
# The last option gives the output file for the combined environment response. The output is saved in 2-bytes integers with 0=voids, 2=walls, 3=filaments and 4=clusters. 
nexus_combine = '{0}/MMF_combine.py node_clean.MMF fila_clean.MMF wall_clean.MMF all_clean.MMF'.format( binDir )




# The following command computes the geometrical direction of the filaments
nexus_direction_1 =  '{0}/MMF_spineDirection fila_clean.MMF fila_spineDirection.MMF --nodeFile node_clean.MMF --radius 1 --convergenceFraction 0.01'.format( binDir )

# The following command computes the geometrical direction of walls (the normal to the plane of the wall)
nexus_direction_2 =  '{0}/MMF_spineDirection wall_clean.MMF wall_spineDirection.MMF --nodeFile node_clean.MMF --filaFile fila_clean.MMF --radius 1 --convergenceFraction 0.01'.format( binDir )

# The following command computes the diameter and linear mass density of filaments
nexus_properties_1 =  '{0}/MMF_spineProperties fila_spineDirection.MMF fila_spineProperties.MMF --densityFile {1} --cleanFile fila_clean.MMF --nodeFile node_clean.MMF --radius 1 --convergenceFraction 0.01'.format( binDir, densityFile )

# The following command computes the thickness and surface mass density of walls
nexus_properties_2 =  '{0}/MMF_spineProperties wall_spineDirection.MMF wall_spineProperties.MMF --densityFile {1} --cleanFile wall_clean.MMF --nodeFile node_clean.MMF --filaFile fila_clean.MMF --radius 1 --convergenceFraction 0.01'.format( binDir, densityFile )




# Writes all the commands in a list. Later on a loop over the list will execute all commands in the specified order.
commands = [
            nexus_response_1,
            nexus_threshold_1,
            nexus_response_2,
            nexus_threshold_2,
            nexus_threshold_3,
            nexus_combine,
            # nexus_direction_1,
            # nexus_direction_2,
            # nexus_properties_1,
            # nexus_properties_2
            ]







#
# No need to modify below this -- this executes the loop commands and checks that each command was run succesfully. It stops if it finds an error.
#
import os
import sys
import numpy as np

noParts = len(commands)
minPart, maxPart = 0, noParts
if len(sys.argv)>=3:
    minPart, maxPart = int(sys.argv[1]), int(sys.argv[2])+1
    if minPart>noParts: sys.exit(1)
    if maxPart>noParts: maxPart=noParts

outputLog = ''
try:
    outputLog = ' >> %s' % logFile
except:
    pass


# run the code for each part
for i in range(minPart,maxPart):    #loop over all parts
    print("\tRunning part %i of %i parts ..." % (i,noParts))
    
    toRun = commands[i] + outputLog
    print("\nRUNNING:\n\t%s\n" % toRun)
    
    status = os.system( toRun )
    if status != 0: # an error took place
        sys.exit(1)
    
    
    

