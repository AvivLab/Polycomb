cd /polycomb
num=2020_233
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.5 0.5 4.0 0.2 2020_233
python3 upload.py polycomb $num results/$num
