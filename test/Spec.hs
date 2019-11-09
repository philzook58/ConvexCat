import LinRel

r20 :: HLinRel IV IV
r20 = resistor 20

main :: IO ()
main = do   putStrLn "Test suite not yet implemented"
            print (r20 == (hid <<< r20))
            print (r20 == r20 <<< hid)
            print (r20 == (hmeet r20 r20))
            print $ resistor 50 == r20 <<< (resistor 30)
            print $ (bridge 10) <<< (bridge 10) == (bridge 5)
            print $ v2h (h2v r20) == r20
            print $ r20 <= htop
            print $ hconverse (hconverse r20) == r20
            print $ (open <<< r20) == open







