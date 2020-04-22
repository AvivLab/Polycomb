cd /polycomb
num=2020_194
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 0.5 6.0 0.1 2020_194
python3 upload.py polycomb $num results/$num
