#Create a softlink for JAseC program.
fn=joinpath(ENV["HOME"], "JAseC")
if isfile(fn)
    run(`unlink $fn`)
end
run(`ln -s $(@__DIR__)/../bin/JAseC $fn`)
