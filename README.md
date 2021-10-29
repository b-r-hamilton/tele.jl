# tele

## To make package
Run in Julia (anywhere)
```
t = Template(; user="b-r-hamilton", dir="~/Code", authors="Brynnydd Hamilton", julia=v"1.6")
```
and then ```t("tele.jl")```
This will create a package template in Code with the name tele. It will also create a GitHub repo with the name "tele". I had to go onto Github website and manually change the name to "tele.jl" to be able to push. Do a force push to overwrite the current data
**Do not have a hyphen in your package name!!!**
