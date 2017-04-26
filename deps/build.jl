using BinDeps
@BinDeps.setup

libhpmpc = library_dependency("libhpmpc")
commit = "cdbaa59518e6a966b01a002c83848c332c47fcdd"

provides(Sources, URI("https://github.com/giaf/hpmpc/archive/$commit.tar.gz"),
         [libhpmpc], unpacked_dir="hpmpc-$commit")

prefix = joinpath(BinDeps.depsdir(libhpmpc), "usr")
srcdir = joinpath(BinDeps.srcdir(libhpmpc), "hpmpc-$commit")

provides(SimpleBuild,
         (@build_steps begin
            GetSources(libhpmpc)
            CreateDirectory(joinpath(prefix, "lib"))
            @build_steps begin
              ChangeDirectory(srcdir)
              `make shared_library USE_BLASFEO=0`
              `mv libhpmpc.$(Libdl.dlext) $prefix/lib/libhpmpc.$(Libdl.dlext)`
            end
          end), [libhpmpc], os=:Unix)

@BinDeps.install Dict(:libhpmpc => :libhpmpc)
