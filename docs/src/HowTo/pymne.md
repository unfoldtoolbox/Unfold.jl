# Loading Data into Unfold

Unfold is generally agnostic to how you load your data. You only require a Matrix (channel x time) or 3D-Array(channel x time x epochs) and an event-dataframe.

## Setup

```@Example main
using Unfold
using UnfoldMakie,CairoMakie

using DataFrames
```

## MNE Demo Dataset

The easiest way to showcase this is to simply use a demo-dataset from MNE.

To extract the metadata, we need `Pandas` which is not installed by default, we therefore need to install it via `CondaPkg`
```@Example main
using CondaPkg
CondaPkg.add("pandas")
CondaPkg.add("mne") # due to a bug in PyMNE https://github.com/beacon-biosignals/PyMNE.jl/issues/38 - we have to add mne additionally
using PyMNE
```

Now we are ready to load some data.

```@Example main
limo_epochs = PyMNE.datasets.limo.load_data(subject=1,path="~/MNE/DATA",update_path=false)
limo_epochs
```

After loading, we can fit an `Unfold` model to it.

First extract the data & convert it to Julia/Unfold requirements



```@Example main
data = pyconvert(Array,limo_epochs.get_data(picks="B11"))
data  = permutedims(data,[2,3,1]) # permute to ch x times x epochs Array format

events = DataFrame(PyTable(limo_epochs.metadata))
rename!(events,2 => :coherence)
events.face = string.(events.face)

times = pyconvert(Vector,limo_epochs.times)
```

Next fit an Unfold Model:

```@Example main
uf = fit(UnfoldModel,[Any=>(@formula(0~face+coherence),times)],events,data)
results = coeftable(uf)
```

```@Example main
plot_results(results)
```

## Read some of your own data

We can make use of all PyMNE importer functions to load the data. Try it for your own data! Get starting with Unfold in no-time!

```@Example main
#eeglabdata = PyMNE.io.read_raw_eeglab("pathToEEGLabSet.set")
```

## Contribute?

Some extra conversions are needed to import the data from PyMNE to Unfold (as shown above). We could try putting these in a wrapper function - do you want to tackle this challenge? Would be a great first contribution to the toolbox :-)
