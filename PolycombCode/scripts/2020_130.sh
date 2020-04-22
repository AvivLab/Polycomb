cd /polycomb
num=2020_130
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.3 2.0 1.0 0.1 2020_130
python3 upload.py polycomb $num results/$num
