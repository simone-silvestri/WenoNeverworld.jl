using WenoNeverworld
using JLD2

dir_path = "/storage4/WenoNeverworldData/"

filename = "weno_thirtytwo_compressed_iteration"
remove_last_character(s) = s[1:end-1]

files = readdir(dir_path)
files = filter((x) -> length(x) >= length(filename), files)
myfiles = filter((x) -> x[1:length(filename)] == filename, files)

myfiles = remove_last_character.(myfiles)
numbers = parse.(Int, filter.(isdigit, myfiles))

@show myfiles
@show numbers
for num in numbers
   new_filename = "weno_thirtytwo_compressed_iteration_new" * string(num) * ".jld2"
   old_filename = "weno_thirtytwo_compressed_iteration" * string(num) * ".jld2"
   
   old_file = jldopen(old_filename)

   jldopen(joinpath(dir_path, new_filename), "w") do f
      for k in keys(old_file)
         if k ∈ ["clock", "resolution"]
            f[k] = old_file[k * "/data"]
         else
            @show k
            f[k * "/data"] = old_file[k * "/data"]
          
         end
	 print(num, " ", k, "\n")
      end
   end
end
