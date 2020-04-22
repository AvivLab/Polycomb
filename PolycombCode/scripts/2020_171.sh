cd /polycomb
num=2020_171
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.5 0.5 6.0 0.15 2020_171
python3 upload.py polycomb $num results/$num
