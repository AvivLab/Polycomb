cd /polycomb
num=2020_175
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.3 0.5 4.0 0.15 2020_175
python3 upload.py polycomb $num results/$num