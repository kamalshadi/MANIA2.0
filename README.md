# To Dos
### Acronyms:
 - __NAR__: Normalized assymetry ratio  
 - __NOS__: Number of streamline  
### 1. Wide Plateau  

Why do we see wide flat region in the __Threshold vs NAR__ plot?
The first assumeption was that in this region the netwrok is stable and not changing by much. But density plot in this region clearly show that network is changing linearly to the threshold. Thus there should be another phenomeno causing this plateau.

> My assumption is that we have alternating network edge change. Because if we add only bidirectional eges clrealy the asymmetry should fall and it should rise in the other scenario, thus we should alternating effect here - At threshold ```T```  a new unidirection edge is added and at ```T+e``` a bidirectional edge and this alternating addition continues in the plateau.

### 2. Envelope of the seed to ROI distance

In analysis of ROI to ROI distance, we have four categories of seed in the seed ROIs.
 - Clearly not connected  
 - Outlier  
 - Envelope seeds  
 - The rest  

The clearly not connected seeds retain a low NOS to the target ROI no matter the distance. The envelope seeds are the ones which in each distance bin have the largest NOS to the target. The outliers are the ones that resides above/below the envelope.

> The question is that can we find this envelope in a robust manner such that it can be parametrized with a smooth mathematical function and if so would this model repeat across ROIs and subjects?

### 3. Distance correction

Depending on the question in 2, we can follow two pathes:
1 - If we have such function we can correct distance bias at the vey begining by using the interplotaed function model and all seed to target ROI connection numbers we feed to MANIA will be already accounted for distance.  

2 - If finding such function is not possible, we talked about two version of the distance correction idea  
 - **The original Idea:** applying the MANIA in distance bins for the subset of ROIs that falls within that bin  
 - **MANIA at seed Level:** Applying MANIA in a __group of seeds__ to taget ROI mode. That is each seed ROI is broken down into groups adn each group has realtively similar distance to the target ROI of choice.  

# CSV File added

The file **L1.csv** has two columns:

1. All distances from seeds in L1 to L4  
2. All connection strengths from seeds in L1 to L4  

File **L4.csv** has the same structure seeding form __L1__


# The list of subjects and their processing status

**L is for Left hemisphere and R for right hemisphere**  

| Subject | Release | Acquisition | Gender | Age | Probtrackx | MANIA | Note |  
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |  
| 188145	|S1200	|Q12|	F |	31-35 |  Running |  |  |  
| 192237	|S1200	|Q11|	F |	31-35 | Done | L Done | re-run for traget ROI 180 |  
| 206323	|S1200	|Q12|	F |	31-35 | Done || re-run for traget ROI 180 |  
| 206929	|S1200	|Q12|	M |	31-35 | Done | L Done | re-run for traget ROI 180 |  
| 208630	|S1200	|Q12|	M |	31-35 | Done || re-run for traget ROI 180 |  
| 227533	|S1200	|Q08|	F |	31-35 | Done | L Done | re-run for traget ROI 180 |  
| 238033	|S1200	|Q09|	M |	31-35 | Done || re-run for traget ROI 180 |  
| 248238	|S1200	|Q13|	F |	31-35 | Done ||re-run for traget ROI 180 |  
| 325129	|S1200	|Q13|	M |	31-35 |  |  |  |  
| 360030	|S1200	|Q13|	F |	31-35 | Done |  |  |  
| 368753	|S1200	|Q12|	F |	31-35 | Done |  |  |  
| 392447	|S1200	|Q07|	M |	31-35 | Done |  |  |  
| 401422	|S1200	|Q08|	F |	31-35 | Done |  |  |   
| 413934	|S1200	|Q09|F | 31-35 | Done |  |  |   
| 421226	|S1200	|Q13|	M |	31-35 | Done |  |  |  
| 453542	|S1200	|Q12|	F |	31-35 | Done |  |  |   
| 454140	|S1200	|Q12|	M |	31-35 | Running |  |  |   
| 463040	|S1200	|Q13|	F |	31-35 | Running |  |  |   
| 468050	|S1200	|Q08|	F |	31-35 | Running |  |  |   
| 481042	|S1200	|Q12|	F |	31-35 | Running |  |  |   
