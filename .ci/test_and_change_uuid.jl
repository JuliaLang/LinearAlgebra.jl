@static if Base.VERSION >= v"1.6"
    using TOML
    using Test
else
    using Pkg: TOML
    using Test
end

# To generate the new UUID, we simply modify the first character of the original UUID
const original_uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
const new_uuid      = "47e2e46d-f89d-539d-b4ee-838fcccc9c8e"

# `@__DIR__` is the `.ci/` folder.
# Therefore, `dirname(@__DIR__)` is the repository root.
const project_filename = joinpath(dirname(@__DIR__), "Project.toml")

@testset "Test that the UUID is unchanged" begin
    project_dict = TOML.parsefile(project_filename)
    @test project_dict["uuid"] == original_uuid
end

write(
    project_filename,
    replace(
        read(project_filename, String),
        r"uuid = .*?\n" => "uuid = \"$(new_uuid)\"\n",
    ),
)
