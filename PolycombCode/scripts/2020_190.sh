cd /polycomb
num=2020_190
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.5 0.5 1.0 0.15 2020_190
python3 upload.py polycomb $num results/$num
