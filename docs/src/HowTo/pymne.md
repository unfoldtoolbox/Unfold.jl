# Loading Data into Unfold

Unfold is generally agnostic to how you load your data. You only require a Matrix (channel x time) or 3D-Array(channel x time x epochs) and an event-dataframe.

## Setup

```@Example main
using Unfold
using UnfoldMakie,CairoMakie
using PyMNE
using DataFrames
```

## MNE Demo Dataset

The easiest way to showcase this is to simply use a demo-dataset from MNE.

```@Example main
limo_epochs = PyMNE.datasets.limo.load_data(subject=1,path="~/MNE/DATA",update_path=false)
limo_epochs
```

Now we can fit a simple `Unfold` model to it.

First extract the data & convert it to Julia/Unfold requirements

```@Example main
data = limo_epochs.get_data(picks="B11")
data  = permutedims(data,[2,3,1]) # get into ch x times x epochs

function convert_pandas(df_pd)
      df= DataFrame()
    for col in df_pd.columns
        df[!, col] = getproperty(df_pd, col).values
    end
    return df
end
events = convert_pandas(limo_epochs.metadata)
rename!(events,2=>:coherence) # negative signs in formulas are not good ;)
events.face = string.(events.face) # ugly names, but fast

```

Next fit an Unfold Model

```@Example main
uf = fit(UnfoldModel,[Any=>(@formula(0~face+coherence),Float64.(limo_epochs.times))],events,data)
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
