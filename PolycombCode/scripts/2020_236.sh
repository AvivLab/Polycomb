cd /polycomb
num=2020_236
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.3 0.5 1.0 0.2 2020_236
python3 upload.py polycomb $num results/$num
