main = do
line <- getLine
let (l:r:_) = map read $ words line
putStrLn $ drop (pred l) $ take r "atcoder"