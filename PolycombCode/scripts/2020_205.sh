cd /polycomb
num=2020_205
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.4 0.5 4.0 0.1 2020_205
python3 upload.py polycomb $num results/$num
