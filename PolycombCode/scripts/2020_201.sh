cd /polycomb
num=2020_201
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.3 0.5 4.0 0.1 2020_201
python3 upload.py polycomb $num results/$num
