cd /polycomb
num=2020_102
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.3 2.0 1.0 0.15 2020_102
python3 upload.py polycomb $num results/$num
