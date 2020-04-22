cd /polycomb
num=2020_57
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.3 1.0 6.0 0.2 2020_57
python3 upload.py polycomb $num results/$num
