cd /polycomb
num=2020_227
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.3 0.5 4.0 0.2 2020_227
python3 upload.py polycomb $num results/$num
