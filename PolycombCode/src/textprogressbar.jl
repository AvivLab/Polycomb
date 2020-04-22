#@debug function textprogressbar(c,strCR="n")
function textprogressbar(c,strCR)
# This function creates a text progress bar. It should be called with a
# STRING argument to initialize and terminate. Otherwise the number correspoding
# to progress in # should be supplied.
# INPUTS:   C   Either: Text string to initialize or terminate
#                       Percentage number to show progress
# OUTPUTS:  N/A
# Example:  Please refer to demo_textprogressbar.m

# Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
# Version: 1.0
# Changes tracker:  29.06.2010  - First version

# Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

    ## Initialization
    #persistent strCR           #   Carriage return pesistent variable

    # Vizualization parameters
    strPercentageLength = 10   #   Length of percentage string (must be >5)
    strDotsMaximum      = 10   #   The total number of dots in a progress bar

    ## Main

    if isempty(strCR) && ~isa(c,String)
        # Progress bar must be initialized with a string
        error("The text progress must be initialized with a string")
    elseif isempty(strCR) && isa(c,String)
        # Progress bar - initialization
        @printf("%s",c)
        strCR = -1
    elseif ~isempty(strCR) && isa(c,String)
        # Progress bar  - termination
        strCR = []
        print("$c\n")
    elseif isa(c, Number)
        # Progress bar - normal progress
        c = floor(c)
        percentageOut = "$c%"
        percentageOut = string(percentageOut,repeat(" ", strPercentageLength-length(percentageOut)-1))
        nDots = convert(Int64,floor(c/100*strDotsMaximum))
        dotOut = string("[", repeat("=",nDots), repeat(" ",strDotsMaximum-nDots),"]")
        strOut = string(percentageOut,dotOut)

        # Print it on the screen
        if strCR == -1
            # Don't do carriage return during first run
            print(strOut)
        else
            # Do it during all the other runs
            print(strCR, strOut)
        end

        # Update carriage return
        strCR = repeat("\b", length(strOut))
        #@bp
    else
        # Any other unexpected input
        error("Unsupported argument type")
    end

end
