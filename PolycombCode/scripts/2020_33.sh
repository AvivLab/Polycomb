cd /polycomb
num=2020_33
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.4 1.0 6.0 0.1 2020_33
python3 upload.py polycomb $num results/$num
