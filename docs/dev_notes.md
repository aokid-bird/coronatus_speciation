# NOTES ON THE SNAKEMAKE
## Ways to pass parameters to R scripts
There are two ways to pass parameters from snakemake to R scripts
### 1: script:
By specifying Rscript in the script: section (Note that the path is relative to the rules somehow, so you need to do "../scripts/"), you can pass any parameters to R and inside R, you can simply do...
```r
param1 <- snakemake@params[['param1']]
param2 <- snakemake@config[['raxml']][['bs']]
param3 <- snakemake@threads
```
### 2. shell & R --args
Use Rscript with --args option in the terminal, like...
```bash
R --vanilla --slave --args {params.param1} {config["raxml"]["bs"]} {output.raxml} < scripts/script.R 
```

## How to use wildcards in general rules that can be applied to any inputs/outputs generated from similar rules
If you want to make a rule general while you want to embed within series of specific rules, then the use of wildcards is useful. For example, if you make a rule "bcf2vcf" that simply defines a function to convert a bcf file to a vcf, and you want to use it in any types of ANGSD runs, then you can do the following.
1. In the rule, for the wildcard, use ```lambda``` in ```input:```, and use ```{{}}``` in ```output:```. This is because only ```input:``` allows you to use a function. In ```output:```, snakemake will understand it based on dependencies from ```input:```.
```python
input:
   bcf = lambda wc: f"results/angsd_{wc.angsd_run}/{output_prefix}/gl.bcf"
output:
    vcf = f"results/angsd_{{angsd_run}}/{output_prefix}/gl.vcf.gz"
```
2. Then you should specify what the wildcard "angsd_run" would indicate in ```rule: all:```. Here, you use expand (the snakemake function) for ```{{angsd_run}}```. Note that the double "{}" is required for escape because it will be first read by python to open the first "{}" (i.e., interpret ```output_prefix``` passed from the python environment), and then it will be passed to expand (snakemake function) where angsd_run is still wrapped by "{}". Because ```angsd_run``` is expanded in the list ```angsd_run = [""]```, you can pass several variables in this, for example ```["raxml", "global"]``` if you need to run this bcf2vcf for multiple sets of ANGSD runs. You can also generalize this by making the part ```angsd_{{angsd_run}}``` to ```{{analysis}}```, where you can specify any analysis name inside ```expand``` 
```python
expand(f"results/angsd_{{angsd_run}}/{output_prefix}/gl.vcf.gz", angsd_run = ["raxml"])
```

## Cautions on conda: section ##
Do not use conda: but include "conda run -n envname" after
#nvname is established by
"conda env create --quiet --file workflow/envs/{env}.yaml"
({env} should include the file name)
at the front-end due to internet inaccessiblility

## Cautions on container: section ##
Do no use container: "docker://url" but include
container: "workflow/containers/{container}.sif"
({caontainer} should include the file name)
after pulling a sif file of the container in the local
at the front-end due to internet inaccessiblility

## Cautions on the use of f-string ##
 ~Before Snakemake 8.2.3, f-string cannot be handled correctly under Python 3.12 or related environments (https://github.com/snakemake/snakemake/issues/2648). This can be resolved by updating Snakemake or avoiding the use of f-string.~ 

NOTES: Somehow, f-string started to work (possibly due to the cache in the .snakemake dir?). And Snakemake >v.8 does not allow the use of --cluster, but requires profiles/pbs/config.yaml with --profile, which did not work correctly due to some conflicts with plugins and PBS. Therefore, it is better to stick with Snakemake v.7. Even if f-string cannot be used, format and "" +var+ "" notation should be used.

In the current snakemake.yaml, the snakemake environment (bioinfo_pipeline) will be created with python 3.10 and Snakemake 7.32.4.

### Usage of f-string with expand
When to use expand with f-string, use {{}} to wrap the variables to expand and use {} for normal wildcards for f-string.
```python
# e.g.
genos = expand(f"results/angsd_intersect/{output_prefix}/{{group}}/gl.geno.gz", group=groups)
```

# NOTES ON GITHUB
## When conflicts occurred between different versions of push requests
When changes were pushed without pulling the previous changes (normally done in the different environment), it will cause errors since there are two diveregent branches of the same name. In this case, one possible solution could be the following
```bash
git pull --rebase origin test
```
This will allow you to pull the different version first, comparing it with your new version, and if there is no conflicts between the changes you've made and the divergent branch, then it will combine the divergent branches.

# NOTES ON APPTAINER/SINGULARITY
If you need softwares that you need to compile by yourself, they are not prepared in Conda channels, and their containers are not available online (DockerHub, BioContainers), then you need to build your own containers. Containers can be built using Apptainer (previously, singularity), which can only be used in Linux environment. If you want to do this in your Mac system, then Virtual Linux Environment should be prepared. 

```lima``` is the best virtual environment to do so.
## IMPORTANT NOTES
Because the architectures are different between Ubuntu (amd x84_64) and Mac Apple Silicon (arm), the sif image built on a normal apptainer Lima virtual environment cannot be used in the Ubuntu amd environment. This causes a problem when you want to execute the built sif image in the PBS cluster. The solution is to build a virtual environment by lima with x86_64 architecture.

https://qiita.com/ikegamitky/items/f40f54e4a7efe1c73daf

I followed this link to create a x84_64 Ubuntu environment, where apptainer was used to build a sif image.


## 1 Set up lima environment
Follow the link above to create the specific environment

4. Check where lima folders are mounted in your Mac machine.
Either of the below will work, which will show you folders mounted to the MacMachine
```bash
mount | grep 9p
mount | grep -i fuse
mount | grep -i virtio
``` 
Then, find the directory with "rw" (read-write), and find it on your Mac. This will be the location where you can share files between your VM and Mac. It is often the case that /Users/Username is read-only, so you can put your .def file in your project directory under your Mac Home, but you can only make .sif in the mounted directory (e.g., /tmp/lima), so you need to move your generated .sfi files by your own.

## 2 Generate a .sif file based on your .def file
1. Generate a sif file
Inside your Linux VM, you build .sif by
```bash
apptainer build /tmp/lima/glactools.sif /path/to/glactools.def
```
where /tmp/lima is your VM mount shared with the host Mac, and /path/to/.def is where you store your def files in your Mac (check if the path is at least read-only from VM).
2. In your Mac terminal, move the sif file to the project dir
```bash
mv /tmp/lima/glactools.sif /path/to/project/glactools.sif
```
/path/to/project/ could be like, gallinago/workflow/containers/glactools/glactools.sif
3. checksum your sif file
```bash
sha256sum /path/to/glactools.sif > /path/to/glactools.sif.sha256
```

## 3. Upload the generated sif and its checksum to the PBS cluster
Send the sif file to the cluster environment
