cd /polycomb
num=2020_93
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.3 2.0 4.0 0.15 2020_93
python3 upload.py polycomb $num results/$num
