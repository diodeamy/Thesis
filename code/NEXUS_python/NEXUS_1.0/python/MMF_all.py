#!/usr/bin/python


import os
import sys
import optparse


def programOptions():
    """ Sets the program options and reads them from the command line. """
    parser = optparse.OptionParser(usage='%prog  scriptFile  [options]')
    parser.add_option('-l','--log', action='store', default=None, dest="logFile", metavar='logFile', help="Choose this option if you would like to send the program messages to a log file. Give the name of the log file.")
    parser.add_option('--all', action='store_true', default=False, dest="allOn", help="Choose this option if to compute the MMF response for all environments (i.e. nodes, filaments and walls).")
    parser.add_option('--node', action='store_true', default=False, dest="nodeOn", help="Choose this option if to compute the MMF response for node environments.")
    parser.add_option('--fila', action='store_true', default=False, dest="filaOn", help="Choose this option if to compute the MMF response for filament environments.")
    parser.add_option('--wall', action='store_true', default=False, dest="wallOn", help="Choose this option if to compute the MMF response for wall environments.")
    parser.add_option('--responseAll', action='store_true', default=False, dest="responseAll", help="Choose this option if to compute the maximum response for all environments during the same call to the maximum response function.")
    
    options, args = parser.parse_args()
    if len(args) <1:
        print "This script needs at least one option giving the name of the script file used to execute the commands."
        sys.exit(1)
    options.scriptFile = args[0]
    if not os.path.isfile(options.scriptFile):
        print "Could not find the script file '%s'. Please check that you gave the correct name for the script file!" % options.scriptFile
        sys.exit(1)
    options.args = args[1:]
    
    if options.allOn: options.nodeOn, options.filaOn, options.wallOn = True, True, True
    
    return options


def run(command):
    print "RUNNING: \n\t" "%s \n" % command
    status = os.system( command )
    if status != 0:
        print "\n~~~~ ERROR ~~~~ while executing the script commands. Check the script output for the reason of the error."
        sys.exit(1)


options = programOptions()
additionalArguments = " ".join(options.args)
if options.logFile is not None:
    additionalArguments += " >> %s" % options.logFile


# compute the response for all environments at once
if options.responseAll:
    tempDict = { 1:"node", 10:"fila", 100:"wall", 11:"node fila", 101:"node wall", 111:"node fila wall" }
    choice = tempDict[ int(options.nodeOn*1+options.filaOn*10+options.wallOn*100) ]
    run( "%s response %s %s" % (options.scriptFile,choice,additionalArguments) )
    if options.nodeOn:
        run( "%s threshold clean node %s" % (options.scriptFile,additionalArguments) )
    if options.filaOn:
        run( "%s threshold clean fila %s" % (options.scriptFile,additionalArguments) )
    if options.wallOn:
        run( "%s threshold clean wall %s" % (options.scriptFile,additionalArguments) )
    sys.exit(0)

# execute the commands for nodes
if options.nodeOn:
    run( "%s response node %s" % (options.scriptFile,additionalArguments) )
    run( "%s threshold clean node %s" % (options.scriptFile,additionalArguments) )

# execute the commands for nodes
if options.filaOn:
    run( "%s response fila %s" % (options.scriptFile,additionalArguments) )
    run( "%s threshold clean fila %s" % (options.scriptFile,additionalArguments) )

# execute the commands for nodes
if options.wallOn:
    run( "%s response wall %s" % (options.scriptFile,additionalArguments) )
    run( "%s threshold clean wall %s" % (options.scriptFile,additionalArguments) )
