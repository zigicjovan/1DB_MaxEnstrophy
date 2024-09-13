%{
check testcase name before execution:
main_1DBEnstrophy_2024_v*(gradient method,desired name for testcase (if duplicate gradient method test))
%}

for testnumber = 2:20    
    for s_w = 1:3
        main_1DBEnstrophy_2024_v4_5('RSWCGJ','4',s_w,testnumber)
        main_1DBEnstrophy_2024_v4_5('RSWCGJ','3',s_w,testnumber)
        main_1DBEnstrophy_2024_v4_5('RSWCGJ','5',s_w,testnumber)
        
        main_1DBEnstrophy_2024_v4_5('RSWG','4',s_w,testnumber)
        main_1DBEnstrophy_2024_v4_5('RSWG','3',s_w,testnumber)
        main_1DBEnstrophy_2024_v4_5('RSWG','5',s_w,testnumber)
        
        main_1DBEnstrophy_2024_v4_5('RSWM','4',s_w,testnumber)
        main_1DBEnstrophy_2024_v4_5('RSWM','3',s_w,testnumber)
        main_1DBEnstrophy_2024_v4_5('RSWM','5',s_w,testnumber)
    end
end