import Pkg; Pkg.activate(".")

using Mangal
using DataFrames
using CSV

total_networks = count(MangalNetwork)
count_per_page = 100
number_of_pages = convert(Int, ceil(total_networks/count_per_page))

LS = DataFrame(fill(Int64, 5), [:id, :S, :L, :P, :H], total_networks)

global cursor = 1
@progress "Paging networks" for page in 1:number_of_pages
    global cursor
    networks_in_page = Mangal.networks("count" => count_per_page, "page" => page-1)
    @progress "Counting items" for current_network in networks_in_page
        S = count(MangalNode, current_network)
        L = count(MangalInteraction, current_network)
        P = count(MangalInteraction, current_network, "type" => "predation")
        H = count(MangalInteraction, current_network, "type" => "herbivory")
        LS[cursor,:] .= (current_network.id, S, L, P, H)
        cursor = cursor + 1
    end
end

CSV.write("ls.csv", LS)
