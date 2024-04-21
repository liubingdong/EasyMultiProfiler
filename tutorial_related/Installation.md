

## EasyMultiProfiler: An Efficient and Convenient R package in Multi-omics Down-Stream Analysis and Visualization
<a href="https://github.com/liubingdong/EasyMultiProfier/blob/main/man/figures/logo.png"><img src="https://github.com/liubingdong/EasyMultiProfier/blob/main/man/figures/logo.png" width=150 align="right" ></a>
![](https://img.shields.io/badge/R%20language->=4.3-brightgreen.svg)
![](https://img.shields.io/badge/Mac%20OSX%20&%20Windows-Available-brightgreen.svg)
![](https://img.shields.io/badge/Release%20version-0.1.0-brightgreen.svg)
[![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/liubingdong/EasyMultiProfier)

The EasyMultiProfiler package aims to offer a user-friendly and efficient multi-omics data analysis tool on the R platform. It facilitates various essential tasks related to microbiome, genome, and metabolite downstream analysis, providing a seamless workflow from start to finish.

### Installation issues

**Easily install**
```R
if (!requireNamespace("pak", quietly=TRUE)) install.packages("pak")
pak::pkg_install("liubingdong/EasyMultiProfiler")
library(EasyMultiProfiler)
```
#### 1. Common error about the timeout
Solution method:
>For some region with unstable network, users could utilize the local mirrors to avoid unexperted errors before installation.
```R
## For china main land users could use this
local({r <- getOption("repos")
r["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
options(repos=r)}
)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("download.file.method"="libcurl")
options("url.method"="libcurl")
```
#### 2. Common error about the <u>pkgbuild::check build tools(debug = TRUE)</u>

Solution method:

>For Windows users may encounter an error <u>"Could not find tools necessary to compile a package"<u> during the installation process. To address this, it's essential to install Rtools beforehand (For R 4.3.x need rtool43, for R 4.4.x need rtool44, [click here ~ 400MB](https://mirrors.tuna.tsinghua.edu.cn/CRAN/)). Afterward, simply restat R and re-try ```pak::pkg_install("liubingdong/EasyMultiProfiler")```.

<img src="Installation_figs/rtool.jpg" alt="rtool" style="zoom:40%;" />

**Completely install** 
```R
# In the step, please type in : 1 2 3 4 5 6 7 
setRepositories(addURLs = c(BioCsoft = "https://bioconductor.org/packages/3.18/bioc",
                  BioCann = "https://bioconductor.org/packages/3.18/data/annotation"))  
options(timeout = 600000000) 
install.packages("remotes") # remotes (>= 2.5.0)
install.packages("BiocManager") # BiocManager (>= 1.30.22)
BiocManager::install("base64enc") # base64enc (>= 0.1.3)
BiocManager::install("WGCNA") # WGCNA (>= 1.72.5)
BiocManager::install("clusterProfiler") # clusterProfiler (>= 4.10.0)
remotes::install_github("liubingdong/EasyMultiProfiler")
library(EasyMultiProfiler)
```
#### 1. Common error about the timeout.

Solution method:
>Because EasyMultiProfiler incluld many necessary data to provide comprehensive tools, the package size (~10MB) may lead to install timeout error for some users in bad network region. Users could reset the Max linktime to make sure the installation success.
```R
options(timeout = 600000000) 
```

#### 2. Common error about the can't find the bioconductor repository by install_github

Solution method:

>Users could set the correct repository handly to make the installation.

```R
setRepositories(addURLs = c(BioCsoft = "https://bioconductor.org/packages/3.18/bioc",
                  BioCann = "https://bioconductor.org/packages/3.18/data/annotation"))  
```

<img src="Installation_figs/setRepositories.jpg" alt="setRepositories" style="zoom:100%;" />

#### 3. Common error about patchwork

Because the patchwork has two version to lead unexpercted conflict in the package enrichplot, users need to degrade the patchwork from 2.4 to 1.2.0.

<img src="Installation_figs/patchwork_error1.jpg" alt="patchwork_error2" style="zoom:100%;" />

Solution method 1:

>During the installation process, be careful that not upgrade patchwork to version 2.4.

Solution method 2:

>If conflict already occuered, users could check the patchwork version and re-install patchwork at the version 1.2.0.

<img src="Installation_figs/patchwork_error2.jpg" alt="patchwork_error2" style="zoom:40%;" />

```R
remotes::install_version("patchwork",version='1.2.0',force = TRUE)
```






