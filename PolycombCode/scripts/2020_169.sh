cd /polycomb
num=2020_169
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.4 0.5 6.0 0.15 2020_169
python3 upload.py polycomb $num results/$num
