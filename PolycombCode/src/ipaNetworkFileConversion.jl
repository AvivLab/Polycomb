# File created by Maryl Lambros on Aug. 26, 2019 to read and convert
# relationship files downloaded from IPA (ingenuity pathway analysis) to use
# to create W matrices for individual that are based on biological reality

using DelimitedFiles

# Function to convert Ipa network relationship files to matrix that can use
# to generate W matrix for code:
function convertIpaNetworkFiles(nameOfFile::String)
    x = joinpath("..","ipaNetworkFiles",nameOfFile)
    data = readdlm(x,'\t')
    data = data[setdiff(1:end, 1:2), :] # remove first two rows that just contain text about file
    data = data[:,setdiff(1:end, 4:5)] # remove last two columns that are empty
    return data
end


# Generate W matrix when choose to use IPA generated gene regulatory network
# relationships file:
function generateWmatUsingIpaNetwork(ipaFileName::String)
    testIpaConnectivity = convertIpaNetworkFiles(ipaFileName)
    allTfs = vcat(intersect(unique(testIpaConnectivity[:,1]),unique(testIpaConnectivity[:,3])),setdiff(unique(testIpaConnectivity[:,1]),unique(testIpaConnectivity[:,3])),setdiff(unique(testIpaConnectivity[:,3]),unique(testIpaConnectivity[:,1])))
    wMatrix = zeros(Float64,length(allTfs),length(allTfs))
    wConnectivityMat = zeros(Float64,length(allTfs),length(allTfs))
    for i = 1:size(testIpaConnectivity)[1]
        tf = testIpaConnectivity[i,1]
        tfIndex = findall(x -> x == tf, allTfs)[1]
        tfRegulee = testIpaConnectivity[i,3]
        tfReguleeIndex = findall(x -> x == tfRegulee, allTfs)[1]
        if testIpaConnectivity[i,2] == "activation"
            wConnectivityMat[tfReguleeIndex,tfIndex] = 1
            wMatVal = randn()
            while wMatVal < 0
                wMatVal = randn()
            end
            wMatrix[tfReguleeIndex,tfIndex] = wMatVal
        elseif testIpaConnectivity[i,2] == "inhibition"
            wConnectivityMat[tfReguleeIndex,tfIndex] = -1
            wMatVal = randn()
            while wMatVal > 0
                wMatVal = randn()
            end
            wMatrix[tfReguleeIndex,tfIndex] = wMatVal
        end
    end
    return wMatrix, wConnectivityMat
end
