// stolen from baz, including this comment vv
// date magic stolen directly from stackoverflow, lol
Date.prototype.isLeapYear = function() {
    const year = this.getFullYear();
    if ((year & 3) != 0) return false;
    return (((year % 100)+100)%100 != 0 || ((year % 400)+400)%400 == 0);
};

Date.prototype.getDOY = function() {
    const dayCount = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    const mn = this.getMonth();
    const dn = this.getDate();
    let dayOfYear = dayCount[mn] + dn;
    if (mn > 1 && this.isLeapYear()) dayOfYear++;
    return dayOfYear;
};

function padDigits(argamStr, digits) {
    let [integral, decimal] = argamStr.split(".");

    // decimal is undefined if there weren't any decimal points
    decimal = (decimal ?? "").padEnd(digits, "0");

    return `${integral}.${decimal}`;
}

function daysBetween(start_y, start_m, start_d, end_y, end_m, end_d) {
    // The number of milliseconds in all UTC days (no DST)
    const oneDay = 1000 * 60 * 60 * 24;

    // A day in UTC always lasts 24 hours (unlike in other time formats)
    const start = Date.UTC(start_y, start_m, start_d);
    const end = Date.UTC(end_y, end_m, end_d);

    // so it's safe to divide by 24 hours
    return Math.round((end-start) / oneDay);
}

function dateToNum(y,m,d){
    return daysBetween(1970,0,1,y,m,d);
}

function numToDate(num){
    new_date = new Date(Date.UTC(1970,0,1+num));
    return [new_date.getUTCFullYear(),new_date.getUTCMonth(),new_date.getUTCDate()];
}

function cardYearFromGregYear(greg_year) {
    // returns the card year that ends in the given greg year
    // A year 0 in the card calendar ends in 1980 CE
    return greg_year-1980;
}

function doesGregYearHaveJoker(greg_year) {
    const arcana = ((cardYearFromGregYear(greg_year)%22)+22)%22;
    if(arcana==6||arcana==11||arcana==17){
        return true;
    }
    return (arcana==0&&((greg_year%20)+20)%20!=0);
}

function cardYearStartDate(card_year) {
    meta_cycle_count = Math.floor(card_year/220);
    
    year_in_meta_cycle = ((card_year%220)+220)%220;
    cycle_count = Math.floor(year_in_meta_cycle/22);
    
    year_in_cycle = ((year_in_meta_cycle%22)+22)%22;
    jokers_this_cycle = 0;
    if(cycle_count>0&&year_in_cycle>0){
        jokers_this_cycle += 1;
    }
    if(year_in_cycle>6){
        jokers_this_cycle += 1;
    }
    if(year_in_cycle>11){
        jokers_this_cycle += 1;
    }
    if(year_in_cycle>17){
        jokers_this_cycle += 1;
    }
    days_since_card_epoch = (
                            meta_cycle_count*(22*52+3+(22*52+4)*9)+
                            Math.max(0,cycle_count-1)+(cycle_count*(22*52+3))+
                            52*year_in_cycle+jokers_this_cycle
                            )*7;
    return dateToNum(1979,8,23)+days_since_card_epoch;
}

function cardDateFromGregDate(greg_year,greg_month,greg_day) {
    card_year_that_ends = cardYearFromGregYear(greg_year);
    year_end_UTC_date = cardYearStartDate(card_year_that_ends+1);
    greg_UTC_date = dateToNum(greg_year,greg_month,greg_day);
    if(greg_UTC_date>=year_end_UTC_date){
        current_card_year = card_year_that_ends+1;
    }else{
        current_card_year = card_year_that_ends;
    }
    days_into_card_year = (greg_UTC_date-cardYearStartDate(current_card_year));
    weeks_into_card_year = Math.floor(days_into_card_year/7);
    seasons_into_card_year = Math.floor(weeks_into_card_year/13);
    week_of_season = ((weeks_into_card_year%13)+13)%13;
    day_of_week = Math.round(((days_into_card_year%7)+7)%7);
    return [current_card_year,seasons_into_card_year,week_of_season,day_of_week];
}

function gregDateFromCardDate(card_year,card_season,card_week,card_day) {
    start_date = cardYearStartDate(card_year);
    return numToDate(start_date+(card_season*13+card_week)*7+card_day);
}

function cardCalenderSVG(card_year){
    days_of_week = ['S','M','T','W','T','F','S'];
    season_names = ['Spades','Diamonds','Clubs','Hearts'];
    season_symbols = ['&spades;','&diams;','&clubs;','&hearts;','&starf;'];
    week_names = ['A','2','3','4','5','6','7','8','9','10','J','Q','K'];
    bg_color = 'FFFFFF';
    black_colors = ['999999','B7B7B7','CCCCCC','D9D9D9'];
    red_colors = ['EB0000','FF3E3E','FF7575','FF8989'];
    major_arcana = ['the Fool',
'the Magician',
'the High Priestess',
'the Empress',
'the Emperor',
'the Hierophant',
'the Lovers',
'the Chariot',
'Strength',
'the Hermit',
'the Wheel of Fortune',
'Justice',
'the Hanged Man',
'Death',
'Temperance',
'the Devil',
'the Tower',
'the Star',
'the Moon',
'the Sun',
'Judgement',
'the World'];
    greg_month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
    color_sets = [black_colors,red_colors];
    s = [24,36,52];

    total_width = ((days_of_week.length+2)*s[1]+s[0])*season_names.length+s[0];
    total_height = s[0]*5+s[1]+s[2]+week_names.length*s[0];

    start_date = numToDate(cardYearStartDate(card_year));
    start_year = start_date[0];
    has_joker = doesGregYearHaveJoker(start_year+1);
    if(has_joker){
        total_height += s[0]*4+s[1];
    }

    svg_out = '<svg width="'+(total_width).toString()+'" height="'+(total_height).toString()+'">\n\t';
    svg_out = svg_out+'<defs>\n\t\t<style type="text/css">@import url("https://fonts.googleapis.com/css?family=Montserrat")\n\t\t</style>\n\t</defs>\n\t';
    svg_out = svg_out+'<rect x="0" y="0" width="'+(total_width).toString()+'" height="'+(total_height).toString()+'" stroke-opacity="0" fill="#'+bg_color+'"/>\n\t';
    for(var i=0;i<season_names.length;i+=1){
        svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])).toString()+'" y="'+(s[0]*2+s[2]).toString()+'" width="'+(s[1]*(2+days_of_week.length)).toString()+'" height="'+(s[1]).toString()+'" stroke-opacity="0" fill="#'+color_sets[i%2][0]+'"/>\n\t';
        svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])).toString()+'" y="'+(s[0]*2+s[2]+s[1]).toString()+'" width="'+(s[1]*(2+days_of_week.length)).toString()+'" height="'+(s[0]*(2+week_names.length)).toString()+'" stroke-opacity="0" fill="#'+color_sets[i%2][1]+'"/>\n\t';
        for(var j=0;j<week_names.length;j+=1){
            svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]).toString()+'" y="'+(s[0]*2+s[2]+s[1]+s[0]*(1+j)).toString()+'" width="'+(s[1]*days_of_week.length).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[i%2][3-(j%2)]+'"/>\n\t';
        }
    }
    if(has_joker){
        for(var i=0;i<days_of_week.length+2;i+=1){
            svg_out = svg_out+'<rect x="'+(s[0]+(season_names.length-1)/2*((days_of_week.length+2)*s[1]+s[0])+s[1]*i).toString()+'" y="'+(s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[0]).toString()+'" width="'+(s[1]).toString()+'" height="'+(s[1]).toString()+'" stroke-opacity="0" fill="#'+color_sets[i%2][0]+'"/>\n\t';
            svg_out = svg_out+'<rect x="'+(s[0]+(season_names.length-1)/2*((days_of_week.length+2)*s[1]+s[0])+s[1]*i).toString()+'" y="'+(s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[0]+s[1]).toString()+'" width="'+(s[1]).toString()+'" height="'+(s[0]*(2+1)).toString()+'" stroke-opacity="0" fill="#'+color_sets[i%2][1]+'"/>\n\t';
        }
        for(var i=1;i<days_of_week.length+1;i+=1){
            for(var j=0;j<1;j+=1){
                svg_out = svg_out+'<rect x="'+(s[0]+(season_names.length-1)/2*((days_of_week.length+2)*s[1]+s[0])+s[1]*i).toString()+'" y="'+(s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[1]+s[0]+s[0]*(1+j)).toString()+'" width="'+(s[1]).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[i%2][3-(j%2)]+'"/>\n\t';
            }

        }
    }
    year_name = 'Deck of '+major_arcana[((card_year%22)+22)%22]+' ('+(start_year).toString()+'-'+(start_year+1).toString()+')';
    svg_out = svg_out+'<text x="'+(total_width/2).toString()+'" y="'+(3+s[0]+s[2]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="36">'+year_name+'</text>'
    for(var i=0;i<season_names.length;i+=1){
        svg_out = svg_out+'<text x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(s[1]*(2+days_of_week.length))/2).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+season_names[i]+' Season</text>'
        svg_out = svg_out+'<text x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(s[1]/2)).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+season_symbols[i]+'</text>'
        svg_out = svg_out+'<text x="'+(s[0]+(i+1)*((days_of_week.length+2)*s[1]+s[0])-s[0]-s[1]+(s[1]/2)).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+season_symbols[i]+'</text>'
        for(var j=0;j<days_of_week.length;j+=1){
            svg_out = svg_out+'<text x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(j+1)*s[1]+s[1]/2).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]+s[0]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="16">'+days_of_week[j]+'</text>'
        }
        for(var j=0;j<week_names.length;j+=1){
            svg_out = svg_out+'<text x="'+(-5+s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="end" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="16">'+week_names[j]+'</text>'
            for(var k=0;k<days_of_week.length;k+=1){
                greg_date = gregDateFromCardDate(card_year,i,j,k);
                date_string = greg_date[2].toString();
                if(greg_date[2]==1){
                    date_string = '1 '+greg_month[greg_date[1]];
                }
                svg_out = svg_out+'<text x="'+(4+s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(k+1)*s[1]).toString()+'" y="'+(-2+s[0]*2+s[2]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="start" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="10">'+date_string+'</text>'
            }
        }
    }
    if(has_joker){
        svg_out = svg_out+'<text x="'+(s[0]+(season_names.length-1)/2*((days_of_week.length+2)*s[1]+s[0])+(s[1]*(2+days_of_week.length))/2).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[0]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">Joker Week</text>'
        svg_out = svg_out+'<text x="'+(s[0]+(season_names.length-1)/2*((days_of_week.length+2)*s[1]+s[0])+(s[1]/2)).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[0]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+season_symbols[4]+'</text>'
        svg_out = svg_out+'<text x="'+(s[0]+(season_names.length-1+2)/2*((days_of_week.length+2)*s[1]+s[0])-s[0]-s[1]+(s[1]/2)).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[0]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+season_symbols[4]+'</text>'
        for(var j=0;j<days_of_week.length;j+=1){
            svg_out = svg_out+'<text x="'+(s[0]+(season_names.length-1)/2*((days_of_week.length+2)*s[1]+s[0])+(j+1)*s[1]+s[1]/2).toString()+'" y="'+(2+s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[0]+s[1]+s[0]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="16">'+days_of_week[j]+'</text>'
        }
        for(var j=0;j<1;j+=1){
            //svg_out = svg_out+'<text x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]).toString()+'" y="'+(s[0]*2+s[2]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="end" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="16">'+week_names[j]+'</text>'
            for(var k=0;k<days_of_week.length;k+=1){
                greg_date = gregDateFromCardDate(card_year,4,j,k);
                date_string = greg_date[2].toString();
                if(greg_date[2]==1){
                    date_string = '1 '+greg_month[greg_date[1]];
                }
                svg_out = svg_out+'<text x="'+(4+s[0]+(season_names.length-1)/2*((days_of_week.length+2)*s[1]+s[0])+(k+1)*s[1]).toString()+'" y="'+(-2+s[0]*2+s[2]+s[1]+s[0]*(2+week_names.length)+s[0]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="start" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="10">'+date_string+'</text>'
            }
        }
    }
    //<rect x="24" y="136" width="324" height="360" stroke-opacity="0" fill="#B7B7B7"/>
    svg_out = svg_out+'</svg>';
    return svg_out;
}

function gregCalenderSVG(greg_year){
    days_of_week = ['S','M','T','W','T','F','S'];
    season_names = ['Spades','Diamonds','Clubs','Hearts'];
    season_symbols = ['&spades;','&diams;','&clubs;','&hearts;','&starf;'];
    week_names = ['A','2','3','4','5','6','7','8','9','10','J','Q','K'];
    bg_color = 'FFFFFF';
    black_colors = ['999999','B7B7B7','CCCCCC','D9D9D9'];
    red_colors = ['EB0000','FF3E3E','FF7575','FF8989'];
    major_arcana = ['the Fool',
'the Magician',
'the High Priestess',
'the Empress',
'the Emperor',
'the Hierophant',
'the Lovers',
'the Chariot',
'Strength',
'the Hermit',
'the Wheel of Fortune',
'Justice',
'the Hanged Man',
'Death',
'Temperance',
'the Devil',
'the Tower',
'the Star',
'the Moon',
'the Sun',
'Judgement',
'the World'];
    greg_month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
    greg_month_full = ['January','February','March','April','May','June','July','August','September','October','November','December'];
    color_sets = [black_colors,red_colors];
    s = [24,36,52];

    start_of_year = Date.UTC(greg_year,0,1);
    end_of_year = Date.UTC(greg_year+1,0,1);


    for(var i=0;i<daysBetween(greg_year,0,1,greg_year+1,0,1);i+=1){
        current_date = Date.UTC(greg_year,0,1+i);
    }

    total_width = ((days_of_week.length+2)*s[1]+s[0])*season_names.length+s[0];
    total_height = s[0]*2+s[2]+3*(9*s[0]+s[1]);

    //start_date = numToDate(cardYearStartDate(card_year));
    start_year = greg_year;

    svg_out = '<svg width="'+(total_width).toString()+'" height="'+(total_height).toString()+'">\n\t';
    svg_out = svg_out+'<defs>\n\t\t<style type="text/css">@import url("https://fonts.googleapis.com/css?family=Montserrat")\n\t\t</style>\n\t</defs>\n\t';
    svg_out = svg_out+'<rect x="0" y="0" width="'+(total_width).toString()+'" height="'+(total_height).toString()+'" stroke-opacity="0" fill="#'+bg_color+'"/>\n\t';
    for(var y=0;y<3;y+=1){
        for(var x=0;x<4;x+=1){
            i = x;
            m = y*4+x;
            month_start = cardDateFromGregDate(greg_year,m,1);
            month_end = cardDateFromGregDate(greg_year,m+1,0);
            console.log(month_start);
            console.log(month_end);
            svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]).toString()+'" width="'+(s[1]*(2+days_of_week.length)).toString()+'" height="'+(s[1]).toString()+'" stroke-opacity="0" fill="#'+color_sets[month_start[1]%2][0]+'"/>\n\t';
            svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]).toString()+'" width="'+(s[1]*(2+days_of_week.length)).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[month_start[1]%2][1]+'"/>\n\t';
            
            greg_start = gregDateFromCardDate(month_end[0],month_end[1],month_end[2],0);
            greg_end = gregDateFromCardDate(month_start[0],month_start[1],month_start[2],0);
            num_weeks = (dateToNum(greg_start[0],greg_start[1],greg_start[2])-dateToNum(greg_end[0],greg_end[1],greg_end[2]))/7+1;
            for(var j=0;j<num_weeks;j+=1){
                push_forward = 0;
                if(j==0){
                    push_forward = month_start[3];
                }
                pull_back = 0;
                if(j==num_weeks-1){
                    pull_back = 6-((((month_end[3]+1)-1)%7)+7)%7;
                }
                week_start = cardDateFromGregDate(greg_year,m,1+7*j);
                week_name_use = week_names[week_start[2]];
                if(week_start[1]==4){
                    for(var k=0;k<days_of_week.length+2;k+=1){
                        svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]*k).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(1+j)).toString()+'" width="'+(s[1]).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[k%2][1]+'"/>\n\t';
                        if(k>push_forward&&k<8-pull_back){
                            svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]*k).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(1+j)).toString()+'" width="'+(s[1]).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[k%2][3-(j%2)]+'"/>\n\t';
                        }
                    }
                }else{
                    svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(1+j)).toString()+'" width="'+(s[1]*(2+days_of_week.length)).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[(month_start[1]+Math.floor((month_start[2]+j)/13))%2][1]+'"/>\n\t';
                    svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]+s[1]*push_forward).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(1+j)).toString()+'" width="'+(s[1]*days_of_week.length-s[1]*(push_forward+pull_back)).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[(month_start[1]+Math.floor((month_start[2]+j)/13))%2][3-(j%2)]+'"/>\n\t';
                }
            }
            svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(num_weeks+1)).toString()+'" width="'+(s[1]*(2+days_of_week.length)).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[month_end[1]%2][1]+'"/>\n\t';
            
        }
    }
    year_name = "Gregorian Calendar Year "+(start_year).toString();
    svg_out = svg_out+'<text x="'+(total_width/2).toString()+'" y="'+(3+s[0]+s[2]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="36">'+year_name+'</text>'
    
    for(var y=0;y<3;y+=1){
        for(var x=0;x<4;x+=1){
            i = x;
            m = y*4+x;
            month_start = cardDateFromGregDate(greg_year,m,1);
            month_end = cardDateFromGregDate(greg_year,m+1,0);
            svg_out = svg_out+'<text x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(s[1]*(2+days_of_week.length))/2).toString()+'" y="'+(2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+greg_month_full[m]+'</text>'
            svg_out = svg_out+'<text x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(s[1]/2)).toString()+'" y="'+(2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+season_symbols[month_start[1]]+'</text>'
            svg_out = svg_out+'<text x="'+(s[0]+(i+1)*((days_of_week.length+2)*s[1]+s[0])-s[0]-s[1]+(s[1]/2)).toString()+'" y="'+(2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="20">'+season_symbols[month_start[1]]+'</text>'
            for(var j=0;j<days_of_week.length;j+=1){
                svg_out = svg_out+'<text x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(j+1)*s[1]+s[1]/2).toString()+'" y="'+(2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]/2).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="middle" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="16">'+days_of_week[j]+'</text>'
            }
            greg_start = gregDateFromCardDate(month_end[0],month_end[1],month_end[2],0);
            greg_end = gregDateFromCardDate(month_start[0],month_start[1],month_start[2],0);
            num_weeks = (dateToNum(greg_start[0],greg_start[1],greg_start[2])-dateToNum(greg_end[0],greg_end[1],greg_end[2]))/7+1;
            console.log(num_weeks);
            //num_weeks = (((month_end[2]-month_start[2]+1)%13)+13)%13;
            for(var j=0;j<num_weeks;j+=1){
                // push_forward = 0;
                // if(j==0){
                //     push_forward = month_start[3];
                // }
                // pull_back = 0;
                // if(j==(((month_end[2]-month_start[2]+1)%13)+13)%13-1){
                //     pull_back = 6-((((month_end[3]+1)-1)%7)+7)%7;
                // }
                // svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(1+j)).toString()+'" width="'+(s[1]*(2+days_of_week.length)).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[(month_start[1]+Math.floor((month_start[2]+j)/13))%2][1]+'"/>\n\t';
                // svg_out = svg_out+'<rect x="'+(s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]+s[1]*push_forward).toString()+'" y="'+(y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(1+j)).toString()+'" width="'+(s[1]*days_of_week.length-s[1]*(push_forward+pull_back)).toString()+'" height="'+(s[0]).toString()+'" stroke-opacity="0" fill="#'+color_sets[(month_start[1]+Math.floor((month_start[2]+j)/13))%2][3-(j%2)]+'"/>\n\t';
                

                week_start = cardDateFromGregDate(greg_year,m,1+7*j);
                week_name_use = week_names[week_start[2]];
                if(week_start[1]==4){
                    week_name_use = season_symbols[4];
                }
                svg_out = svg_out+'<text x="'+(-5+s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]).toString()+'" y="'+(2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="end" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="16">'+week_name_use+'</text>'
                for(var k=0;k<days_of_week.length;k+=1){
                    greg_date = gregDateFromCardDate(month_start[0],month_start[1],month_start[2]+j,k);
                    date_string = greg_date[2].toString();
                    if(greg_date[1]==m){
                        svg_out = svg_out+'<text x="'+(4+s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(k+1)*s[1]).toString()+'" y="'+(-2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="start" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="10">'+date_string+'</text>'
                    }
                }
            }
            // for(var j=0;j<week_names.length;j+=1){
            //     svg_out = svg_out+'<text x="'+(-5+s[0]+i*((days_of_week.length+2)*s[1]+s[0])+s[1]).toString()+'" y="'+(2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="end" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="16">'+week_names[(((month_start[2]+j)%13)+13)%13]+'</text>'
            //     for(var k=0;k<days_of_week.length;k+=1){
            //         greg_date = gregDateFromCardDate(greg_year,i,j,k);
            //         date_string = greg_date[2].toString();
            //         if(greg_date[2]==1){
            //             date_string = '1 '+greg_month[greg_date[1]];
            //         }
            //         svg_out = svg_out+'<text x="'+(4+s[0]+i*((days_of_week.length+2)*s[1]+s[0])+(k+1)*s[1]).toString()+'" y="'+(-2+y*(9*s[0]+s[1])+s[0]*2+s[2]+s[1]+s[0]*(j+1+0.5)).toString()+'" style="font-family:Montserrat" dominant-baseline="middle" text-anchor="start" fill="#000000" fill-opacity="1" font-weight="bold" font-style="bold" font-size="10">'+date_string+'</text>'
            //     }
            // }
        }
    }
    //<rect x="24" y="136" width="324" height="360" stroke-opacity="0" fill="#B7B7B7"/>
    svg_out = svg_out+'</svg>';
    return svg_out;
}

function set_initial_time() {
    const date = new Date();

    const year = date.getFullYear();
    const month = date.getMonth();
    const day = date.getDate();
    const weekday = date.getDay();

    const days_since_card_epoch = daysBetween(2021,8,22,year,month,day); // card epoch starts 2021-09-22
    const time_since_card_epoch = days_since_card_epoch-0.803219697; // except actually it's technically 19:16:38 UTC on that day so account for that

    const toS = (x) => x.toString();
    const toN = (x) => x.toString();

    const gregDateString = toS(year) + '-' + ('0'+(month+1).toString()).slice(-2) + '-' + ('0'+day.toString()).slice(-2);
    document.getElementById("date_input").value = gregDateString;
}
function update_time(date_str) {
    const date = new Date(date_str);
    console.log(date);
    console.log(date.getUTCDate());

    const year = date.getUTCFullYear();
    const month = date.getUTCMonth();
    const day = date.getUTCDate();
    const weekday = date.getUTCDay();

    const days_since_card_epoch = daysBetween(2021,8,22,year,month,day); // card epoch starts 2021-09-22
    const time_since_card_epoch = days_since_card_epoch-0.803219697; // except actually it's technically 19:16:38 UTC on that day so account for that

    const toS = (x) => x.toString();
    const toN = (x) => x.toString();

    const gregDateString = toS(year) + '-' + ('0'+(month+1).toString()).slice(-2) + '-' + ('0'+day.toString()).slice(-2);
    card_date = cardDateFromGregDate(year,month,day);

    days_of_week = ['Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'];
    season_names = ['Spades','Diamonds','Clubs','Hearts'];
    season_symbols = ['&spades;','&diams;','&clubs;','&hearts;','&starf;'];
    week_symbols = ['A','2','3','4','5','6','7','8','9','10','J','Q','K'];
    week_names = ['Ace','Two','Three','Four','Five','Six','Seven','Eight','Nine','Ten','Jack','Queen','King'];
    bg_color = 'FFFFFF';
    black_colors = ['999999','B7B7B7','CCCCCC','D9D9D9'];
    red_colors = ['EB0000','FF3E3E','FF7575','FF8989'];
    major_arcana = ['the Fool',
'the Magician',
'the High Priestess',
'the Empress',
'the Emperor',
'the Hierophant',
'the Lovers',
'the Chariot',
'Strength',
'the Hermit',
'the Wheel of Fortune',
'Justice',
'the Hanged Man',
'Death',
'Temperance',
'the Devil',
'the Tower',
'the Star',
'the Moon',
'the Sun',
'Judgement',
'the World'];
    major_arcana_symbols = ['0','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI'];
    planet_symbols = ['&#9737;','&#9790;','&male;','&#9791;','&#9795;','&female;','&#9796;'];

    
    let cardDateFullString = days_of_week[card_date[3]] + ' in the ' + week_names[card_date[2]] + ' of ' + season_names[card_date[1]] + ' from the deck of ' + major_arcana[((card_date[0]%22)+22)%22];
    if(card_date[1]==4){
        cardDateFullString = days_of_week[card_date[3]] + ' in the Joker from the deck of ' + major_arcana[((card_date[0]%22)+22)%22];
    }
    
    let cardDateString = planet_symbols[card_date[3]] + '-' + week_symbols[card_date[2]] + season_symbols[card_date[1]] + '-' + major_arcana_symbols[((card_date[0]%22)+22)%22];
    if(card_date[1]==4){
        cardDateString = planet_symbols[card_date[3]] + '-JOKER-' + major_arcana_symbols[((card_date[0]%22)+22)%22];
    }

    //const year = date.getFullYear();
    const dayOfYear = date.getDOY(); // Math.ceil((date - new Date(year,0,1)) / 86400000);

    const dateStringDigits = toS(year) + '-' +
                         toS(dayOfYear);

    const dateStringName = toN(year) + ' ' +
                       toN(dayOfYear);

    const unix = +(date.getTime() / 1000).toFixed(2);

    const unixStringDigits = padDigits(toS(unix), 2);
    const unixStringName = toN(Math.floor(unix));

    document.getElementById('gregdate').innerHTML = gregDateString;
    document.getElementById('fullcarddate').innerHTML = cardDateFullString;
    document.getElementById('carddate').innerHTML = cardDateString;

}
//setInterval(update_time, 100000);

set_initial_time();
update_time(document.getElementById("date_input").value);
// svg_string = cardCalenderSVG(44);
// var doc = new DOMParser().parseFromString(svg_string, 'application/xml');
// var el = document.getElementById("svg_flag");
// el.innerHTML = svg_string;

function edit_calendar(year){
    document.getElementById("year_after").innerText = (parseInt(year)+1).toString();
    svg_string = cardCalenderSVG(cardYearFromGregYear(parseInt(year)+1));
    var doc = new DOMParser().parseFromString(svg_string, 'application/xml');
    var el = document.getElementById("card_cal");
    el.innerHTML = svg_string;
    document.getElementById("year_after").innerText = (parseInt(year)+1).toString();
    svg_string = gregCalenderSVG(parseInt(year));
    var doc = new DOMParser().parseFromString(svg_string, 'application/xml');
    var el = document.getElementById("greg_cal");
    el.innerHTML = svg_string;
}

edit_calendar(2023);