# input SF file and the name of the assembly 
SF_file=$1
ASM_Name=$2


# split into each superfamily and merge monomers within 

## splitting superfamilies into subgroups based on ages here https://docs.google.com/spreadsheets/d/1lTTWgoKzFxOIhl0BW1jLhyo3cvAKZKuKg0soHaAlzOE/edit?usp=sharing

## first the newest - dhor SF01 and SF02 

grep -e J3 -e J4 -e J5 -e J6 $SF_file >  SF01.bed
bedtools merge -d 350 -i  SF01.bed | awk '($3-$2) >= 450' >  SF01.merged.bed

grep -e D3 -e D4 -e D5 -e D6 -e D7 -e D8 -e D9 $SF_file >  SF02.bed
bedtools merge -d 350 -i  SF02.bed | awk '($3-$2) >= 450' >  SF02.merged.bed


## now the possibly HOR - can be called as active active SF1, SF2, SF3
grep -e J1 -e J2 $SF_file > SF1.bed
bedtools merge -d 350 -i  SF1.bed | awk '($3-$2) >= 450'  >  SF1.merged.bed

grep  -e D1 -e D2 -e FD $SF_file >  SF2.bed
bedtools merge -d 350 -i  SF2.bed | awk '($3-$2) >= 450' >  SF2.merged.bed

grep -e W1 -e W2 -e W3 -e W4 -e W5 $SF_file >  SF3.bed
bedtools merge -d 350 -i  SF3.bed | awk '($3-$2) >= 450' >  SF3.merged.bed

## now the older - also active in primates SF4, SF5, SF6 
grep Ga $SF_file >  SF4.bed
bedtools merge -d 350 -i  SF4.bed | awk '($3-$2) >= 450' >  SF4.merged.bed

# because this one is a dimer we want to have both of the monomers present - so we merge and overlap 
grep -e R1 -e R2 $SF_file > SF5.bed
grep R1 $SF_file >  SF5R1.bed
grep R2 $SF_file >  SF5R2.bed
bedtools merge -d 600 -i  SF5R1.bed >  SF5R1.merged.bed
bedtools merge -d 600 -i  SF5R2.bed >  SF5R2.merged.bed



bedtools intersect -wo -b SF5R1.merged.bed -a SF5R2.merged.bed  \
    | awk '{print $1 "\t" $2 "\t" $3; print $4 "\t" $5 "\t" $6}'  > SF5.temp.bed

# arrays that are pure R2 are also ok so we include all of those entries 
awk '($3-$2) >= 450' SF5R2.merged.bed >> SF5.temp.bed 

bedtools sort -i SF5.temp.bed  | bedtools merge -d 600 -i stdin | awk '($3-$2) >= 450' > SF5.merged.bed

grep Ha $SF_file >  SF6.bed
bedtools merge -d 350 -i  SF6.bed | awk '($3-$2) >= 450' >  SF6.merged.bed

## now ancient SF7, SF8, SF9, SF10, SF11 
grep Ka $SF_file >  SF7.bed
bedtools merge -d 350 -i  SF7.bed | awk '($3-$2) >= 450' >  SF7.merged.bed


# because this is a dimer we need both present 
grep -e Na -e Oa $SF_file > SF8.bed
grep Oa $SF_file >  SF8Oa.bed
grep Na $SF_file >  SF8Na.bed
bedtools merge -d 600 -i  SF8Oa.bed >  SF8Oa.merged.bed
bedtools merge -d 600 -i  SF8Na.bed >  SF8Na.merged.bed

bedtools intersect -wo -a SF8Oa.merged.bed -b SF8Na.merged.bed \
    | awk '{print $1 "\t" $2 "\t" $3; print $4 "\t" $5 "\t" $6}' \
    | bedtools sort -i stdin | bedtools merge -d 350 -i stdin | awk '($3-$2) >= 450' > SF8.merged.bed

grep Ca $SF_file >  SF9.bed
bedtools merge -d 350 -i  SF9.bed | awk '($3-$2) >= 450' >  SF9.merged.bed

grep Ba $SF_file >  SF10.bed
bedtools merge -d 350 -i  SF10.bed | awk '($3-$2) >= 450' >  SF10.merged.bed

grep Ja $SF_file >  SF11.bed
bedtools merge -d 350 -i  SF11.bed | awk '($3-$2) >= 450' >  SF11.merged.bed


## more ancient SF12, SF13, SF14, SF15, SF16
grep Aa $SF_file >  SF12.bed
bedtools merge -d 350 -i  SF12.bed | awk '($3-$2) >= 450' >  SF12.merged.bed

grep Ia $SF_file >  SF13.bed
bedtools merge -d 350 -i  SF13.bed | awk '($3-$2) >= 450' >  SF13.merged.bed

grep La $SF_file >  SF14.bed
bedtools merge -d 350 -i  SF14.bed | awk '($3-$2) >= 450' >  SF14.merged.bed

grep Fa $SF_file >  SF15.bed
bedtools merge -d 350 -i  SF15.bed | awk '($3-$2) >= 450' >  SF15.merged.bed

grep Ea $SF_file >  SF16.bed
bedtools merge -d 350 -i  SF16.bed | awk '($3-$2) >= 450' >  SF16.merged.bed

## prehistoric SF17 and SF18 
grep Qa $SF_file >  SF17.bed
bedtools merge -d 350 -i  SF17.bed | awk '($3-$2) >= 450' >  SF17.merged.bed


grep -e Pa -e Ta $SF_file >  SF18.bed
bedtools merge -d 350 -i  SF18.bed | awk '($3-$2) >= 450' >  SF18.merged.bed

# Now for each I want to evaluate the coverage and then sort into high and low confidence files 

## unfortunately I need to treat each group of superfamilies seperately so no big loops 
## I'm going to use the ages to structure how I break ties - so newer SFs overpower older SFs 
echo > high_conf.bed
echo > low_conf.bed

## now the possibly HOR - can be called as active active SF1, SF2, SF3
bedtools merge -d 9000 -i  SF1.merged.bed \
    | bedtools coverage -a stdin -b SF1.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "hor(SF1)", "0", ".", $2, $3, "255,102,0", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.6) print >> "high_conf.bed"; else if ($10 >= 0.3 && $10 <= 0.6) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF2.merged.bed \
    | bedtools coverage -a stdin -b SF2.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "hor(SF2)", "0", ".", $2, $3, "255,102,0", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.6) print >> "high_conf.bed"; else if ($10 >= 0.3 && $10 <= 0.6) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF3.merged.bed \
    | bedtools coverage -a stdin -b SF3.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "hor(SF3)", "0", ".", $2, $3, "255,102,0", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.6) print >> "high_conf.bed"; else if ($10 >= 0.3 && $10 <= 0.6) print >> "low_conf.bed"}' 


## now the older - also active in primates SF4, SF5, SF6 

bedtools merge -d 9000 -i  SF4.merged.bed \
    | bedtools coverage -a stdin -b SF4.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon/hor(SF4)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.8) print >> "high_conf.bed"; else if ($10 >= 0.5 && $10 <= 0.8) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF5.merged.bed \
    | bedtools coverage -a stdin -b SF5.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon/hor(SF5)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.8) print >> "high_conf.bed"; else if ($10 >= 0.5 && $10 <= 0.8) print >> "low_conf.bed"}' 

## now we find the active array 

bedtools sort -i high_conf.bed | awk '($3-$2) >= 1500' > high_conf.sorted.bed
bedtools sort -i low_conf.bed | awk '($3-$2) >= 1500' > low_conf.sorted.bed


bedtools subtract -a low_conf.sorted.bed -b high_conf.sorted.bed > arrays.bed
cat high_conf.sorted.bed >> arrays.bed
bedtools sort -i arrays.bed > arrays.sorted.bed 

## Use this set  to identify active 
chromosomes=( $(cut -f1 arrays.sorted.bed | uniq) )

#incase of rerun 
touch array.active.bed
rm array.active.bed

for chrom in ${chromosomes[@]}; do 
    grep -w $chrom arrays.sorted.bed > tmp.bed \
    && awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, ($3-$2)}' tmp.bed | sort -nk 10,10 -r > tmp.sorted.bed \
    && head -n 1 tmp.sorted.bed | awk '{gsub("mon/", "");print}' | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "activeSF_"$4, $5, $6, $7, $8, "153,0,0"}' >> array.active.bed 
done


## process the rest of the Superfamilies that are not active 

bedtools merge -d 9000 -i  SF6.merged.bed \
    | bedtools coverage -a stdin -b SF6.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon/hor(SF6)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.8) print >> "high_conf.bed"; else if ($10 >= 0.5 && $10 <= 0.8) print >> "low_conf.bed"}' 

## now ancient SF7, SF8, SF9, SF10, SF11 
bedtools merge -d 9000 -i  SF7.merged.bed \
    | bedtools coverage -a stdin -b SF7.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF7)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.8) print >> "high_conf.bed"; else if ($10 >= 0.6 && $10 <= 0.8) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF8.merged.bed \
    | bedtools coverage -a stdin -b SF8.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF8)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.8) print >> "high_conf.bed"; else if ($10 >= 0.6 && $10 <= 0.8) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF9.merged.bed \
    | bedtools coverage -a stdin -b SF9.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF9)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.8) print >> "high_conf.bed"; else if ($10 >= 0.6 && $10 <= 0.8) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF10.merged.bed \
    | bedtools coverage -a stdin -b SF10.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF10)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.8) print >> "high_conf.bed"; else if ($10 >= 0.6 && $10 <= 0.8) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF11.merged.bed \
    | bedtools coverage -a stdin -b SF11.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF11)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.6 && $10 <= 0.8) print >> "low_conf.bed"}' 

## more ancient SF12, SF13, SF14, SF15, SF16
bedtools merge -d 9000 -i  SF12.merged.bed \
    | bedtools coverage -a stdin -b SF12.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF12)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.7 && $10 <= 0.9) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF13.merged.bed \
    | bedtools coverage -a stdin -b SF13.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF13)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.7 && $10 <= 0.9) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF14.merged.bed \
    | bedtools coverage -a stdin -b SF14.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF14)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.7 && $10 <= 0.9) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF15.merged.bed \
    | bedtools coverage -a stdin -b SF15.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF15)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.7 && $10 <= 0.9) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF16.merged.bed \
    | bedtools coverage -a stdin -b SF16.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF16)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.7 && $10 <= 0.9) print >> "low_conf.bed"}' 


## prehistoric SF17 and SF18 
bedtools merge -d 9000 -i  SF17.merged.bed \
    | bedtools coverage -a stdin -b SF17.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF17)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.8 && $10 <= 0.9) print >> "low_conf.bed"}' 

bedtools merge -d 9000 -i  SF18.merged.bed \
    | bedtools coverage -a stdin -b SF18.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "mon(SF18)", "0", ".", $2, $3, "255,204,153", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.9) print >> "high_conf.bed"; else if ($10 >= 0.8 && $10 <= 0.9) print >> "low_conf.bed"}' 


## now the newest - dhor SF01 and SF02 - these will have the lowest threshold to be considered high conf 
bedtools merge -d 9000 -i  SF01.merged.bed \
    | bedtools coverage -a stdin -b SF01.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "dhor(SF01)", "0", ".", $2, $3, "255,146,0", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.3) print >> "high_conf.bed"}' 

bedtools merge -d 9000 -i  SF02.merged.bed \
    | bedtools coverage -a stdin -b SF02.bed \
    | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "dhor(SF02)", "0", ".", $2, $3, "255,146,0", $7}' \
    | awk 'BEGIN {OFS="\t"} {if ($10 > 0.3) print >> "high_conf.bed"}' 

# Sort high and low confidence for all of the 
bedtools sort -i high_conf.bed | awk '($3-$2) >= 1500' > high_conf.sorted.bed
bedtools sort -i low_conf.bed | awk '($3-$2) >= 1500' > low_conf.sorted.bed

# First we can finalize the detailed version 
bedtools subtract -a high_conf.sorted.bed -b array.active.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' > arrays.bed
bedtools subtract -a low_conf.sorted.bed -b high_conf.sorted.bed | bedtools subtract -a stdin -b array.active.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' >> arrays.bed
cat array.active.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' >> arrays.bed 
bedtools sort -i arrays.bed > arrays.sorted.bed

#to finalize we just add the rest of the monomeric 

## first just get all monomers 
bedtools merge -d 9000 -i $SF_file >  mon_merged.bed


# remove any that we have already annotated and label the rest Mon
bedtools subtract -a  mon_merged.bed -b  arrays.sorted.bed \
    | awk '($3-$2) >= 350' \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "mon", "0", ".", $2, $3, "255,204,153"}' > summary.bed 

cat arrays.sorted.bed >> summary.bed 

# Finalize the detailed file 
echo 'track name="Detailed_'$ASM_Name'_Alpha_Summary" visibility=2 itemRgb="On"' > ${ASM_Name}_alphaSummary_detailed.sorted.bed
bedtools sort -i summary.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' >> ${ASM_Name}_alphaSummary_detailed.sorted.bed
rm summary.bed

#For the simple summary we are only going to look at relatively large arrays 


## filter to remove anything that is smaller than ~ 3 HOR repeats 6000 
## first in the high confidence
cat high_conf.sorted.bed | awk '($3-$2) >= 6000' | grep hor > high_conf.filtered.bed

## now low confidence - needs to be a very large array
cat low_conf.sorted.bed | awk '($3-$2) >= 10000' | grep hor > low_conf.filtered.bed

# now merge again for the simple summary 
bedtools subtract -a high_conf.filtered.bed -b array.active.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' > arrays.bed
bedtools subtract -a low_conf.filtered.bed -b high_conf.filtered.bed | bedtools subtract -a stdin -b array.active.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' >> arrays.bed
cat array.active.bed | awk  'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' >> arrays.bed 
bedtools sort -i arrays.bed > arrays.sorted.bed

# remove any that we have already annotated and label the rest Mon
bedtools subtract -a  mon_merged.bed -b  arrays.sorted.bed | bedtools merge -d 9000  \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "mon", "0", ".", $2, $3, "255,204,153"}' \
    | awk '($3-$2) >= 350' >>  summary.bed 

cat arrays.sorted.bed >> summary.bed 

## finalize the simple format
echo 'track name="Simple_'$ASM_Name'_Alpha_Summary" visibility=2 itemRgb="On"' > ${ASM_Name}_alphaSummary_simple.sorted.bed
bedtools sort -i summary.bed | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $2, $3, $9}' >> ${ASM_Name}_alphaSummary_simple.sorted.bed


## cleanup
rm SF*bed
rm high_conf.*
rm low_conf.*
rm array*bed