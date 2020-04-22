cd /polycomb
num=2020_99
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.5 2.0 4.0 0.15 2020_99
python3 upload.py polycomb $num results/$num
