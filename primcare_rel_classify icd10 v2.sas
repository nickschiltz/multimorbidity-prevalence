

libname cond "C:\Users\nks8\Box\Data\MEPS\conditions";
libname full "C:\Users\nks8\Box\Data\MEPS\consolidated";


*********************************************************************;
* Bring in 2016-2019 condition file
*********************************************************************;
data cond1;
	set cond.cond2016 (in=a) cond.cond2017 (in=b)
		cond.cond2018 (in=c) cond.cond2019 (in=d);
	if a then YEAR = 2016;
	if b then YEAR = 2017;
	if c then YEAR = 2018;
	if d then YEAR = 2019;
run;

*flag conditions based on Nicholson et al, 2017 (J comorbidity);
data cond2;
	set cond1;
	if ICD10CDX in ("I10" "I11" "I12" "I13" "I14" "I15") then crc_HYPERTEN = 1;
	if ICD10CDX in ("F33" "F40" "F41") then crc_DEPR_ANXIETY = 1;
	if substr(ICD10CDX,1,2) in ("M4" "M7") then crc_MUSCSKELETAL = 1;
	if ICD10CDX in ("M50" "M51" "M52" "M53" "M54" "M60" "M61" "M62" "M63"
				"M65" "M66" "M67" "M68" "M69") then crc_MUSCSKELETAL = 1;
	if ICD10CDX in ("M05" "M13" "M15" "M16" "M17" "M18" "M19") then crc_ARTHRITIS = 1;
	if ICD10CDX in ("M81") then crc_OSTEOPOROSIS = 1;
	if ICD10CDX in ("J40" "J41" "J42" "J43" "J44" "J45" "J46") then crc_ASTHCOPD = 1;
	if ICD10CDX in ("I20" "I25" "I48") then crc_CARDIO=1;
	if substr(ICD10CDX,1,2) in ("I7") then crc_CARDIO = 1;
	if ICD10CDX in ("I05" "I06" "I07" "I08" "I09" "I34" "I35" "I36" "I37" "I38" "I39"
					"I42" "I43" "I50") then crc_HEARTFAIL = 1;
	if ICD10CDX in ("G45" "I62") then crc_STROKE = 1;
	if ICD10CDX in ("K21" "K25" "K29" ) then crc_STOMACH = 1;
	if ICD10CDX in ("K50" "K51" "K52" "K57" "K58") then crc_COLON = 1;
	if ICD10CDX in ("K70" "K71" "K72" "K73" "K74" "K75" "K76" "K77" ) then crc_HEPATITIS = 1;
	if ICD10CDX in ("E10" "E11" "E12" "E13" "E14" ) then crc_DIABETES = 1;
	if ICD10CDX in ("E00" "E01" "E02" "E03" "E04" "E05" "E06" "E07") then crc_THYROID = 1;
	if substr(ICD10CDX,1,1) in ("C") then crc_CANCER = 1;
	if ICD10CDX in ("N18" "N19" ) then crc_KIDNEY = 1;
	if ICD10CDX in ("N03" "N11" "N18" "N20" "N21" "N22" "N23" "N25" "N26" "N27" "N28" "N29" "N50" "N51") then crc_URINARY = 1;
	if substr(ICD10CDX,1,2) in ("N3" "N4") then crc_URINARY = 1;
	if ICD10CDX in ("F00" "F01" "F02" "F03" ) then crc_DEMENTIA = 1;
	if ICD10CDX in ("E78") then crc_HYPERLIP = 1;
	if ICD10CDX in ("E66") then crc_OBESITY = 1;
*treatment flag;
	if sum(of hhnum ipnum obnum opnum rxnum ernum) ge 1 then treated = 1; else treated = 0;
run;

*frequency each condition appears;
proc freq data=cond2;
	table crc_: treated;
run;

*===============================================
* All self-reproted conditions version
*   	for treated prevalence skip farther down
*===============================================;

proc sort data=cond2; by pandupersid year; run;

proc means data=cond2 noprint;
	by pandupersid year;
	var crc_: ;
	output out = cond3 (drop = _TYPE_ _FREQ_)
	max = ;
run;

*frequency each condition appears - by person (may double count if in 2 year panel);
proc freq data=cond3;
	table crc_: ; 
run;

 proc contents data=full.fullyr2019;
run; 
                                                                       

*==========================================================
* Full year consolidated - person level data
*=========================================================;

*bring in full year consolidated person files;
data pers1;
	set full.fullyr2016 (in=a) full.fullyr2017 (in=b)
		full.fullyr2018 (in=c) full.fullyr2019 (in=d);
	if b then year=2017;
	if c then year=2018;
	if d then year=2019;
	keep year pandupersid varpsu varstr perwt: sex
		age16x age17x age18x age19x
		racethx raceax racebx racewx hispanx inscope 
		faminc16-faminc19 unins16-unins19 inscov16-inscov19
		totexp16-totexp19 rthlth31 rthlth42 rthlth53;
run;
proc contents data=pers1; run;

*clean data and drop anyone out of scope, under age 18, or with 0 person-weight;
*inflation-adjustment guidelines: https://meps.ahrq.gov/about_meps/Price_Index.shtml;
	*Option 2, Guideline 1 used: PCE-Health (Table 3, column 1). ;
data pers2;
	set pers1;
	if year = 2016 then do;
		PERSWT = perwt16f;
		AGE = age16x;
		FAMINC = faminc16;
		INSCOV = inscov16;
		UNINS = abs(unins16 - 2);
		TOTEXP = totexp16;
		TOTEXP2019 = totexp16*(110.675/105.449); *inflation adjusted to 2018 dollars;
	end;
	if year = 2017 then do;
		PERSWT = perwt17f;
		AGE = age17x;
		FAMINC = faminc17;
		INSCOV = inscov17;
		UNINS = abs(unins17 - 2);
		TOTEXP = totexp17;
		TOTEXP2019 = totexp17*(110.675/107.250);
	end;
	if year = 2018 then do;
		PERSWT = perwt18f;
		AGE = age18x;
		FAMINC = faminc18;
		INSCOV = inscov18;
		UNINS = abs(unins18 - 2);
		TOTEXP = totexp18;
		TOTEXP2019 = totexp18*(110.675/109.127);
	end;
	if year = 2019 then do;
		PERSWT = perwt19f;
		AGE = age19x;
		FAMINC = faminc19;
		INSCOV = inscov19;
		UNINS = abs(unins19 - 2);
		TOTEXP = totexp19;
		TOTEXP2019 = totexp19;
	end;

	AVGPERSWT = PERSWT/4;
	FEMALE = SEX - 1;

	*create perceived health and poor perceived health variables;
	if rthlth42 ge 1 then SRHEALTH = rthlth42;
		else if rthlth53 ge 1 then SRHEALTH = rthlth53;
		else if rthlth31 ge 1 then SRHEALTH = rthlth31;
		else SRHEALTH = .;
	if SRHEALTH in (4,5) then POORHLTH = 1;
		else if SRHEALTH in (1,2,3) then POORHLTH = 0;
		else POORHLTH = .;
	if age lt 18 then age_delete = 1;
	if inscope ne 1 then ins_delete = 1;
	if perswt le 0 then wt_delete = 1;
	keep year pandupersid varpsu varstr perswt avgperswt female age
		racethx faminc inscov totexp totexp2019 unins srhealth poorhlth 
		age_delete ins_delete wt_delete;
run;

proc freq data=pers2;
	table age_delete; *drop 32246, keep = 93262;
proc freq data=pers2;
	where age_delete ne 1; *drop 3314, keep = 89948;
	table wt_delete;  
proc freq data=pers2;
	where age_delete ne 1 and wt_delete ne 1;
	table ins_delete; *none - all inscope=0 have perswt=0;
run;
*okay to subset because all stratum and psu in the subset;
data pers3;
	set pers2;
	if age_delete = 1 or wt_delete=1 then delete;
	drop age_delete ins_delete wt_delete;
run;

*some data checks;
proc contents data=pers3; run;
proc print data=pers3 (obs=20); run;
proc freq data=pers3;
	table srhealth poorhlth unins;
run;
proc means data=pers3;
	var totexp totexp2019 perswt;
run;

*=================================================================
* Merge full year with conditions file
*================================================================;

proc sql;
	create table cohort1 as
	select *
	from pers3 as a left join cond3 as b
	on (a.pandupersid=b.pandupersid) and (a.year=b.year);
quit;

proc print data=cohort1 (obs=20);
run;

*clean condition data - change . to 0;
data cohort2;
	set cohort1;
	array cond[*] crc_: ;
	do i = 1 to dim(cond);
		if cond[i] = . then cond[i] = 0;
	end;
	drop i;
run;
proc print data=cohort2 (firstobs = 5001 obs=5020); run;

*check for missing -  just n=1 missing for self-rated health;
proc mi data=cohort2 nimpute=0;
run;

*check frequencies of crc for final cohort;
proc freq data=cohort2;
	table crc_: ;
run;

*export to CSV for association rule mining in R;
proc export data=cohort2
   outfile='C:\Users\nks8\Box\General_Reference\M\MEPS KL2 project\assoc_rules\chrcond_2016_2019_v1.csv'
   dbms=csv
   replace;
run;


*=======================================================================
* Treated prevalence version
*======================================================================;

*make one record per person of just treated prevalence;
*----------------------------------------------------;
data cond4t;
	set cond2;
	where treated = 1;
run;

proc sort data=cond4t; by pandupersid year; run;

proc means data=cond4t noprint;
	by pandupersid year;
	var crc_: ;
	output out = cond5t (drop = _TYPE_ _FREQ_)
	max = ;
run;

*merge treated prevalence with person details;
*----------------------------------------------;
proc sql;
	create table cohort1t as
	select *
	from pers3 as a left join cond5t as b
	on (a.pandupersid=b.pandupersid) and (a.year=b.year);;
quit;

proc print data=cohort1t (obs=20);
run;

*clean condition data - change . to 0;
data cohort2t;
	set cohort1t;
	array cond[*] crc_: ;
	do i = 1 to dim(cond);
		if cond[i] = . then cond[i] = 0;
	end;
	drop i;
run;
proc print data=cohort2t (firstobs = 5001 obs=5020); run;

*check for missing -  just n=4 are missing for self-rated health;
proc mi data=cohort2t nimpute=0;
run;

*check frequencies of crc for final cohort;
proc freq data=cohort2t;
	table crc_: ;
run;

*export to CSV for association rule mining in R;
proc export data=cohort2t
   outfile='C:\Users\nks8\Box\General_Reference\M\MEPS KL2 project\assoc_rules\tx_chrcond_2016_2019_v1.csv'
   dbms=csv
   replace;
run;

ods pdf body = "C:\Users\nks8\Box\General_Reference\M\MEPS KL2 project\assoc_rules\treatedprev1619_descriptive.pdf";
*Table 1 - population characteristics;
proc surveymeans data=cohort2t sum mean std median q1 q3;
	var totexp2019;
	cluster VARPSU;
	stratum VARSTR;
	weight AVGPERSWT;
run;

data analysis1;
	set cohort2t;
	*top quartile exp variable;
	MEDEXPQ4 = 0;
	if totexp2019 gt 5813.484 then MEDEXPQ4 = 1;
	*multimorbidity categories;
	NUMCC = sum(of crc_:);
	if NUMCC = 0 then CC0 = 1; else CC0 = 0;
	if NUMCC ge 1 then CC1 = 1; else CC1 = 0;
	if NUMCC ge 2 then CC2 = 1; else CC2 = 0;
	if NUMCC ge 3 then CC3 = 1; else CC3 = 0;
	if NUMCC ge 4 then CC4 = 1; else CC4 = 0;
	if NUMCC ge 5 then CC5 = 1; else CC5 = 0;
	*dummy flag for total count vars;
	DUMMY = 1;
	*age group;
	if 18 le age le 39 then AGEGRP = 1;
	else if 40 le age le 64 then AGEGRP = 2;
	else if age ge 65 then AGEGRP = 3;
run;

proc surveymeans data=analysis1 sum mean std median q1 q3;
	var totexp2019;
	class CC0;
	cluster VARPSU;
	stratum VARSTR;
	weight AVGPERSWT;
run;
 
proc surveyfreq data=analysis1;
	table DUMMY AGEGRP FEMALE RACETHX POORHLTH MEDEXPQ4 CC: / clwt cl;
	cluster VARPSU;
	stratum VARSTR;
	weight AVGPERSWT;
run;

proc contents data=cohort2t;
run;

ods pdf close;






********END PROGRAM*********;




proc sort data=cond3; by dupersid ccscat; run;
*make one record per unique ccscat_system per person;
proc means data=cond3 noprint;
	by dupersid ccscat;
	var chron_adult chron_ped;
	output out=temp7 (drop=_TYPE_ _FREQ_)
	max = chron_adult chron_ped;
	format ccscat FLGDXCCS.;
run;
ods pdf body = "D:\nick_documents\Box Sync\Personal\!current_projects\MEPS KL2 project\sas_output\hwang_ccs_2011-15.pdf";
proc freq data=temp7;
	where chron_adult=1;
	table ccscat;
run;
proc freq data=temp7 order = freq;
	where chron_adult=1;
	table ccscat;
run;
ods pdf close;

proc freq data=cond3;
	where chron_adult=1 and ccscat=657;
	table icd9codx;
run;






	*make rankable categories - 
		*start with OASH (https://www.cdc.gov/pcd/issues/2013/12_0239.htm);
	if ccscat in (98, 99) then rcc_HYPERTN=1; else rcc_HYPERTN=0; *hypertension;
	if ccscat = 108 then rcc_CHF=1; else rcc_CHF=0; *congestive heart failure;
	if ccscat in (100, 101) then rcc_CAD =1; else rcc_CAD=0; *coronary artery disease;
	if ccscat in (105, 106) then rcc_ARYTHM = 1; else rcc_ARYTHM=0; *cardiac arrythmia;
	if ccscat = 53 then rcc_HYPERLIP=1; else rcc_HYPERLIP = 0; *hyperlipedmia;
	if 109 le ccscat le 112 then rcc_STROKE = 1; else rcc_STROKE = 0; *stroke;
	if ccscat in (202,203) then rcc_ARTHR = 1; else rcc_ARTHR = 0; *arthritis;
	if ccscat = 128 then rcc_ASTHMA = 1; else rcc_ASTHMA = 0; *asthma;
	if ccscat = 655 then rcc_AUTISM =1; else rcc_AUTISM = 0; *Autism and others;
	if 11 le ccscat le 43 then rcc_CANCER=1; else rcc_CANCER=0; *cancer;
	if ccscat=108 then rcc_CKD = 1; else rcc_CKD = 0; *chronic kidney disease;
	if ccscat=127 then rcc_COPD = 1; else rcc_COPD = 0; *chr. obstruct pulm disease;
	if ccscat=653 then rcc_DEM = 1; else rcc_DEM = 0; *dementia and Alzh;
	if ccscat=657 then rcc_MOOD=1; else rcc_MOOD = 0; *depression/mood disorders;
	if ccscat in(49,50) then rcc_DIAB=1; else rcc_DIAB=0; *diabetes;

run;

proc freq data=cond3 order=freq;
	table ccscat;
	format ccscat FLGDXCCS.;
run;

proc freq data=cond3 order=freq;
	where chron_adult=1;
	table ccscat;
	format ccscat FDCCSPDX.;
run;
proc freq data=cond3;
	where chron_adult=1;
	table ccscat;
	format ccscat FLGDXCCS.;
run;


proc print data=cond.cond1996 (obs=10);
run;





proc format;
Value FDCCSPDX
           . = ' .: No diagnosis code'
          .A = '.A: Invalid diagnosis'
          .C = '.C: Inconsistent diagnosis'
      1-  10 = ' 1: Infectious and Parasitic DX'
     11-  47 = ' 2: Neoplasms'
     48-  58 = ' 3: Endocr, Nutri, Metab, Immun DX'
     59-  64 = ' 4: Dx of Blood, Blood-Forming Organs'
     65-  75 = ' 5: Mental Health Disorders'
     76-  95 = ' 6: Dx of Nervous System, Sense Organs'
     96- 121 = ' 7: Dx of Circulatory System'
    122- 134 = ' 8: Dx of Respiratory System'
    135- 155 = ' 9: Dx of Digestive System'
    156- 175 = '10: Dx of Genitourinary System'
    176- 196 = '11: Complic Preg, Birth, Puerperium'
    197- 200 = '12: Dx of Skin and Subcutaneous Tissue'
    201- 212 = '13: Dx of Musculoskel, Connective Tissue'
    213- 217 = '14: Congenital Anomalies'
    218- 224 = '15: Perinatal Conditions'
    225- 244 = '16: Injury and Poisoning'
    245- 259 = '17: Other Conditions'
         260 = '18: E codes'
   2601-2621 = '18: E codes'
   651-670	 = ' 5: Mental Health Disorders'
   ;

Value FLGDXCCS
      . = '   .: Missing'
     .A = '   A: Invalid diagnosis'
     .C = '   C: Inconsistent'
     .Z = '   Z: Overall'
      1 = '   1: Tuberculosis'
      2 = '   2: Septicemia (except in labor)'
      3 = '   3: Bacterial infection; unspecified site'
      4 = '   4: Mycoses'
      5 = '   5: HIV infection'
      6 = '   6: Hepatitis'
      7 = '   7: Viral infection'
      8 = '   8: Other infections; including parasitic'
      9 = '   9: Sexually transmitted infections (not HIV or hepatitis)'
     10 = '  10: Immunizations and screening for infectious disease'
     11 = '  11: Cancer of head and neck'
     12 = '  12: Cancer of esophagus'
     13 = '  13: Cancer of stomach'
     14 = '  14: Cancer of colon'
     15 = '  15: Cancer of rectum and anus'
     16 = '  16: Cancer of liver and intrahepatic bile duct'
     17 = '  17: Cancer of pancreas'
     18 = '  18: Cancer of other GI organs; peritoneum'
     19 = '  19: Cancer of bronchus; lung'
     20 = '  20: Cancer; other respiratory and intrathoracic'
     21 = '  21: Cancer of bone and connective tissue'
     22 = '  22: Melanomas of skin'
     23 = '  23: Other non-epithelial cancer of skin'
     24 = '  24: Cancer of breast'
     25 = '  25: Cancer of uterus'
     26 = '  26: Cancer of cervix'
     27 = '  27: Cancer of ovary'
     28 = '  28: Cancer of other female genital organs'
     29 = '  29: Cancer of prostate'
     30 = '  30: Cancer of testis'
     31 = '  31: Cancer of other male genital organs'
     32 = '  32: Cancer of bladder'
     33 = '  33: Cancer of kidney and renal pelvis'
     34 = '  34: Cancer of other urinary organs'
     35 = '  35: Cancer of brain and nervous system'
     36 = '  36: Cancer of thyroid'
     37 = '  37: Hodgkin`s disease'
     38 = '  38: Non-Hodgkin`s lymphoma'
     39 = '  39: Leukemias'
     40 = '  40: Multiple myeloma'
     41 = '  41: Cancer; other and unspecified primary'
     42 = '  42: Secondary malignancies'
     43 = '  43: Malignant neoplasm without specification of site'
     44 = '  44: Neoplasms of unspecified nature or uncertain behavior'
     45 = '  45: Maintenance chemotherapy; radiotherapy'
     46 = '  46: Benign neoplasm of uterus'
     47 = '  47: Other and unspecified benign neoplasm'
     48 = '  48: Thyroid disorders'
     49 = '  49: Diabetes mellitus without complication'
     50 = '  50: Diabetes mellitus with complications'
     51 = '  51: Other endocrine disorders'
     52 = '  52: Nutritional deficiencies'
     53 = '  53: Disorders of lipid metabolism'
     54 = '  54: Gout and other crystal arthropathies'
     55 = '  55: Fluid and electrolyte disorders'
     56 = '  56: Cystic fibrosis'
     57 = '  57: Immunity disorders'
     58 = '  58: Other nutritional; endocrine; and metabolic disorders'
     59 = '  59: Deficiency and other anemia'
     60 = '  60: Acute posthemorrhagic anemia'
     61 = '  61: Sickle cell anemia'
     62 = '  62: Coagulation and hemorrhagic disorders'
     63 = '  63: Diseases of white blood cells'
     64 = '  64: Other hematologic conditions'
     65 = '  65: Mental retardation'
     66 = '  66: Alcohol-related mental disorders'
     67 = '  67: Substance-related mental disorders'
     68 = '  68: Senility and organic mental disorders'
     69 = '  69: Affective disorders'
     70 = '  70: Schizophrenia and related disorders'
     71 = '  71: Other psychoses'
     72 = '  72: Anxiety; somatoform; dissociative; and personality disorders'
     73 = '  73: Preadult disorders'
     74 = '  74: Other mental conditions'
     75 = '  75: Personal history of mental disorder; mental and behavioral problems; observation and screening for mental condition'
     76 = '  76: Meningitis (except that caused by tuberculosis or sexually transmitted disease)'
     77 = '  77: Encephalitis (except that caused by tuberculosis or sexually transmitted disease)'
     78 = '  78: Other CNS infection and poliomyelitis'
     79 = '  79: Parkinson`s disease'
     80 = '  80: Multiple sclerosis'
     81 = '  81: Other hereditary and degenerative nervous system conditions'
     82 = '  82: Paralysis'
     83 = '  83: Epilepsy; convulsions'
     84 = '  84: Headache; including migraine'
     85 = '  85: Coma; stupor; and brain damage'
     86 = '  86: Cataract'
     87 = '  87: Retinal detachments; defects; vascular occlusion; and retinopathy'
     88 = '  88: Glaucoma'
     89 = '  89: Blindness and vision defects'
     90 = '  90: Inflammation; infection of eye (except that caused by tuberculosis or sexually transmitteddisease)'
     91 = '  91: Other eye disorders'
     92 = '  92: Otitis media and related conditions'
     93 = '  93: Conditions associated with dizziness or vertigo'
     94 = '  94: Other ear and sense organ disorders'
     95 = '  95: Other nervous system disorders'
     96 = '  96: Heart valve disorders'
     97 = '  97: Peri-; endo-; and myocarditis; cardiomyopathy (except that caused by tuberculosis or sexually transmitted disease)'
     98 = '  98: Essential hypertension'
     99 = '  99: Hypertension with complications and secondary hypertension'
    100 = ' 100: Acute myocardial infarction'
    101 = ' 101: Coronary atherosclerosis and other heart disease'
    102 = ' 102: Nonspecific chest pain'
    103 = ' 103: Pulmonary heart disease'
    104 = ' 104: Other and ill-defined heart disease'
    105 = ' 105: Conduction disorders'
    106 = ' 106: Cardiac dysrhythmias'
    107 = ' 107: Cardiac arrest and ventricular fibrillation'
    108 = ' 108: Congestive heart failure; nonhypertensive'
    109 = ' 109: Acute cerebrovascular disease'
    110 = ' 110: Occlusion or stenosis of precerebral arteries'
    111 = ' 111: Other and ill-defined cerebrovascular disease'
    112 = ' 112: Transient cerebral ischemia'
    113 = ' 113: Late effects of cerebrovascular disease'
    114 = ' 114: Peripheral and visceral atherosclerosis'
    115 = ' 115: Aortic; peripheral; and visceral artery aneurysms'
    116 = ' 116: Aortic and peripheral arterial embolism or thrombosis'
    117 = ' 117: Other circulatory disease'
    118 = ' 118: Phlebitis; thrombophlebitis and thromboembolism'
    119 = ' 119: Varicose veins of lower extremity'
    120 = ' 120: Hemorrhoids'
    121 = ' 121: Other diseases of veins and lymphatics'
    122 = ' 122: Pneumonia (except that caused by tuberculosis or sexually transmitted disease)'
    123 = ' 123: Influenza'
    124 = ' 124: Acute and chronic tonsillitis'
    125 = ' 125: Acute bronchitis'
    126 = ' 126: Other upper respiratory infections'
    127 = ' 127: Chronic obstructive pulmonary disease and bronchiectasis'
    128 = ' 128: Asthma'
    129 = ' 129: Aspiration pneumonitis; food/vomitus'
    130 = ' 130: Pleurisy; pneumothorax; pulmonary collapse'
    131 = ' 131: Respiratory failure; insufficiency; arrest (adult)'
    132 = ' 132: Lung disease due to external agents'
    133 = ' 133: Other lower respiratory disease'
    134 = ' 134: Other upper respiratory disease'
    135 = ' 135: Intestinal infection'
    136 = ' 136: Disorders of teeth and jaw'
    137 = ' 137: Diseases of mouth; excluding dental'
    138 = ' 138: Esophageal disorders'
    139 = ' 139: Gastroduodenal ulcer (except hemorrhage)'
    140 = ' 140: Gastritis and duodenitis'
    141 = ' 141: Other disorders of stomach and duodenum'
    142 = ' 142: Appendicitis and other appendiceal conditions'
    143 = ' 143: Abdominal hernia'
    144 = ' 144: Regional enteritis and ulcerative colitis'
    145 = ' 145: Intestinal obstruction without hernia'
    146 = ' 146: Diverticulosis and diverticulitis'
    147 = ' 147: Anal and rectal conditions'
    148 = ' 148: Peritonitis and intestinal abscess'
    149 = ' 149: Biliary tract disease'
    150 = ' 150: Liver disease; alcohol-related'
    151 = ' 151: Other liver diseases'
    152 = ' 152: Pancreatic disorders (not diabetes)'
    153 = ' 153: Gastrointestinal hemorrhage'
    154 = ' 154: Noninfectious gastroenteritis'
    155 = ' 155: Other gastrointestinal disorders'
    156 = ' 156: Nephritis; nephrosis; renal sclerosis'
    157 = ' 157: Acute and unspecified renal failure'
    158 = ' 158: Chronic kidney disease'
    159 = ' 159: Urinary tract infections'
    160 = ' 160: Calculus of urinary tract'
    161 = ' 161: Other diseases of kidney and ureters'
    162 = ' 162: Other diseases of bladder and urethra'
    163 = ' 163: Genitourinary symptoms and ill-defined conditions'
    164 = ' 164: Hyperplasia of prostate'
    165 = ' 165: Inflammatory conditions of male genital organs'
    166 = ' 166: Other male genital disorders'
    167 = ' 167: Nonmalignant breast conditions'
    168 = ' 168: Inflammatory diseases of female pelvic organs'
    169 = ' 169: Endometriosis'
    170 = ' 170: Prolapse of female genital organs'
    171 = ' 171: Menstrual disorders'
    172 = ' 172: Ovarian cyst'
    173 = ' 173: Menopausal disorders'
    174 = ' 174: Female infertility'
    175 = ' 175: Other female genital disorders'
    176 = ' 176: Contraceptive and procreative management'
    177 = ' 177: Spontaneous abortion'
    178 = ' 178: Induced abortion'
    179 = ' 179: Postabortion complications'
    180 = ' 180: Ectopic pregnancy'
    181 = ' 181: Other complications of pregnancy'
    182 = ' 182: Hemorrhage during pregnancy; abruptio placenta; placenta previa'
    183 = ' 183: Hypertension complicating pregnancy; childbirth and the puerperium'
    184 = ' 184: Early or threatened labor'
    185 = ' 185: Prolonged pregnancy'
    186 = ' 186: Diabetes or abnormal glucose tolerance complicating pregnancy; childbirth; or the puerperium'
    187 = ' 187: Malposition; malpresentation'
    188 = ' 188: Fetopelvic disproportion; obstruction'
    189 = ' 189: Previous C-section'
    190 = ' 190: Fetal distress and abnormal forces of labor'
    191 = ' 191: Polyhydramnios and other problems of amniotic cavity'
    192 = ' 192: Umbilical cord complication'
    193 = ' 193: OB-related trauma to perineum and vulva'
    194 = ' 194: Forceps delivery'
    195 = ' 195: Other complications of birth; puerperium affecting management of mother'
    196 = ' 196: Normal pregnancy and/or delivery'
    197 = ' 197: Skin and subcutaneous tissue infections'
    198 = ' 198: Other inflammatory condition of skin'
    199 = ' 199: Chronic ulcer of skin'
    200 = ' 200: Other skin disorders'
    201 = ' 201: Infective arthritis and osteomyelitis (except that caused by tuberculosis or sexually transmitted disease)'
    202 = ' 202: Rheumatoid arthritis and related disease'
    203 = ' 203: Osteoarthritis'
    204 = ' 204: Other non-traumatic joint disorders'
    205 = ' 205: Spondylosis; intervertebral disc disorders; other back problems'
    206 = ' 206: Osteoporosis'
    207 = ' 207: Pathological fracture'
    208 = ' 208: Acquired foot deformities'
    209 = ' 209: Other acquired deformities'
    210 = ' 210: Systemic lupus erythematosus and connective tissue disorders'
    211 = ' 211: Other connective tissue disease'
    212 = ' 212: Other bone disease and musculoskeletal deformities'
    213 = ' 213: Cardiac and circulatory congenital anomalies'
    214 = ' 214: Digestive congenital anomalies'
    215 = ' 215: Genitourinary congenital anomalies'
    216 = ' 216: Nervous system congenital anomalies'
    217 = ' 217: Other congenital anomalies'
    218 = ' 218: Liveborn'
    219 = ' 219: Short gestation; low birth weight; and fetal growth retardation'
    220 = ' 220: Intrauterine hypoxia and birth asphyxia'
    221 = ' 221: Respiratory distress syndrome'
    222 = ' 222: Hemolytic jaundice and perinatal jaundice'
    223 = ' 223: Birth trauma'
    224 = ' 224: Other perinatal conditions'
    225 = ' 225: Joint disorders and dislocations; trauma-related'
    226 = ' 226: Fracture of neck of femur (hip)'
    227 = ' 227: Spinal cord injury'
    228 = ' 228: Skull and face fractures'
    229 = ' 229: Fracture of upper limb'
    230 = ' 230: Fracture of lower limb'
    231 = ' 231: Other fractures'
    232 = ' 232: Sprains and strains'
    233 = ' 233: Intracranial injury'
    234 = ' 234: Crushing injury or internal injury'
    235 = ' 235: Open wounds of head; neck; and trunk'
    236 = ' 236: Open wounds of extremities'
    237 = ' 237: Complication of device; implant or graft'
    238 = ' 238: Complications of surgical procedures or medical care'
    239 = ' 239: Superficial injury; contusion'
    240 = ' 240: Burns'
    241 = ' 241: Poisoning by psychotropic agents'
    242 = ' 242: Poisoning by other medications and drugs'
    243 = ' 243: Poisoning by nonmedicinal substances'
    244 = ' 244: Other injuries and conditions due to external causes'
    245 = ' 245: Syncope'
    246 = ' 246: Fever of unknown origin'
    247 = ' 247: Lymphadenitis'
    248 = ' 248: Gangrene'
    249 = ' 249: Shock'
    250 = ' 250: Nausea and vomiting'
    251 = ' 251: Abdominal pain'
    252 = ' 252: Malaise and fatigue'
    253 = ' 253: Allergic reactions'
    254 = ' 254: Rehabilitation care; fitting of prostheses; and adjustment of devices'
    255 = ' 255: Administrative/social admission'
    256 = ' 256: Medical examination/evaluation'
    257 = ' 257: Other aftercare'
    258 = ' 258: Other screening for suspected conditions (not mental disorders or infectious disease)'
    259 = ' 259: Residual codes; unclassified'
    260 = ' 260: E Codes: All (external causes of injury and poisoning)'
    650 = ' 650: Adjustment disorders'
    651 = ' 651: Anxiety disorders'
    652 = ' 652: Attention-deficit, conduct, and disruptive behavior disorders'
    653 = ' 653: Delirium, dementia, and amnestic and other cognitive disorders'
    654 = ' 654: Developmental disorders'
    655 = ' 655: Disorders usually diagnosed in infancy, childhood, or adolescence'
    656 = ' 656: Impulse control disorders, NEC'
    657 = ' 657: Mood disorders'
    658 = ' 658: Personality disorders'
    659 = ' 659: Schizophrenia and other psychotic disorders'
    660 = ' 660: Alcohol-related disorders'
    661 = ' 661: Substance-related disorders'
    662 = ' 662: Suicide and intentional self-inflicted injury'
    663 = ' 663: Screening and history of mental health and substance abuse codes'
    670 = ' 670: Miscellaneous disorders'
   2601 = '2601: E Codes: Cut/pierceb'
   2602 = '2602: E Codes: Drowning/submersion'
   2603 = '2603: E Codes: Fall'
   2604 = '2604: E Codes: Fire/burn'
   2605 = '2605: E Codes: Firearm'
   2606 = '2606: E Codes: Machinery'
   2607 = '2607: E Codes: Motor vehicle traffic (MVT)'
   2608 = '2608: E Codes: Pedal cyclist; not MVT'
   2609 = '2609: E Codes: Pedestrian; not MVT'
   2610 = '2610: E Codes: Transport; not MVT'
   2611 = '2611: E Codes: Natural/environment'
   2612 = '2612: E Codes: Overexertion'
   2613 = '2613: E Codes: Poisoning'
   2614 = '2614: E Codes: Struck by; against'
   2615 = '2615: E Codes: Suffocation'
   2616 = '2616: E Codes: Adverse effects of medical care'
   2617 = '2617: E Codes: Adverse effects of medical drugs'
   2618 = '2618: E Codes: Other specified and classifiable'
   2619 = '2619: E Codes: Other specified; NEC'
   2620 = '2620: E Codes: Unspecified'
   2621 = '2621: E Codes: Place of occurrence';
run;
