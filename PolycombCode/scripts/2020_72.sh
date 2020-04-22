cd /polycomb
num=2020_72
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.5 1.0 4.0 0.2 2020_72
python3 upload.py polycomb $num results/$num
