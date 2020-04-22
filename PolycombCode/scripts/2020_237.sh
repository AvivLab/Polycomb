cd /polycomb
num=2020_237
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.3 0.5 1.0 0.2 2020_237
python3 upload.py polycomb $num results/$num
