cd /polycomb
num=2020_58
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.3 1.0 6.0 0.2 2020_58
python3 upload.py polycomb $num results/$num
