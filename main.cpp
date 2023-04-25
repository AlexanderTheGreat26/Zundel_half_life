#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>


const int T = 100; // Frame period.


typedef std::tuple<int, int, int, int, int> cation; // <Atom ID, atom type, bonds count, Connected atom ID , Connected atom ID>

typedef std::vector<std::vector<std::string>> cations;


cations cations_strings (const std::string & file_name, const std::string & pattern);

void file_creation (const std::string & file_name, cations & stuff, const int & step);

cations Z_remove (const std::string & Z_file_name, const std::string & stuff_file_name, const int & frame_period);


int main() {
    // We need find different Zundel cations. Their sign is string " 2 2 " between atoms numbers.
    cations Z = std::move(cations_strings("bonds.only", " 2 2 "));
    file_creation("Zundels.only", Z, T);
    // Then we need H3O-cations, but in default file we can find parts of Z-cations with pattern " 3 3 ".
    cations stuff = std::move(cations_strings("bonds.only", " 3 3 "));
    file_creation("stuff", stuff, T);
    // Then need remove parts of this Z-cations from corresponding frames and find uniq H3O-cations.
    // For this purpose we will read the stuff-file and use type cation.
    cations H3O = std::move(Z_remove("Zundels.only", "stuff", T));
    file_creation("H3O.only", H3O, T);
    return 0;
}


void file_creation (const std::string & file_name, cations & stuff, const int & step) {
        std::ofstream fout;
        fout.open(file_name, std::ios::trunc);
        for (int i = 0; i < stuff.size()-1; ++i) {
            fout << "#\n# Timestep\t" << i*step << '\n';
            for (const auto & line: stuff[i])
                fout << line << '\n';
        }
        fout.close();
}



bool contain (const std::string & word, const std::string & line) {
    return line.find(word) != std::string::npos;
}


cations cations_strings (const std::string & file_name, const std::string & pattern) {
    cations result;
    std::string line;
    int i = 0;
    std::ifstream fin(file_name);
    if (!fin.is_open()) throw std::runtime_error("Error opening file.");
    while (!fin.eof())
        while (getline(fin, line)) {
            bool timestep_line = contain("Timestep", line);
            if (!(timestep_line || contain("#", line))) {

                if (contain(pattern, line))
                    result[i].emplace_back(line);

            } else if (timestep_line && !result.empty()) {
                ++i;
            } else result.resize(result.size() + 1);
        }
    return result;
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


std::vector<std::string> frame_read (const std::string & file_name, int timestep) {
    std::string line;
    std::vector<std::string> frame;
    bool current_frame = false;

    std::ifstream fin(file_name);
    if (!fin.is_open()) throw std::runtime_error("Error opening " + file_name + "!");
    while (!fin.eof()) {
        while (getline(fin, line)) {
            bool current_step_line = contain("Timestep\t" + toString(timestep), line);

            if (current_step_line) {
                current_frame = true;
                continue;
            }

            if (contain("#\n", line) && current_frame) break;

            if (current_frame) frame.emplace_back(line);

        }
    }
    return frame;
}


cation string2cation (std::string & line) {
    std::istringstream iss(line);
    int atom_ID, atom_type, bond_count, first_bonded_atom_ID, second_bonded_atom_ID;
    iss >> atom_ID >> atom_type >> bond_count >> first_bonded_atom_ID >> second_bonded_atom_ID;
    return std::make_tuple(atom_ID, atom_type, bond_count, first_bonded_atom_ID, second_bonded_atom_ID);
}


void cleaning_stuff (std::vector<cation> & Zundels, std::vector<cation> & stuff) {
    for (auto & Zundel : Zundels)
        for (int j = 0; j < stuff.size(); ++j)
            if (std::get<3>(Zundel) == std::get<0>(stuff[j]) ||
                std::get<4>(Zundel) == std::get<0>(stuff[j]))
                stuff.erase(stuff.begin()+j);
    std::cout << "Here!" << std::endl;
}


template<typename T, size_t... Is>
std::string tuple2string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple2string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple2string_impl(t, std::make_index_sequence<size>{});
}


cations Z_remove (const std::string & Z_file_name, const std::string & stuff_file_name, const int & frame_period) {
    cations clear_H3O;

    std::vector<cation> Z_cations_tplvec, stuff_cations_tplvec;
    int i = 0;

    do {
        clear_H3O.resize(clear_H3O.size()+1);
        Z_cations_tplvec.clear();
        stuff_cations_tplvec.clear();
        std::vector<std::string> Zundels = frame_read(Z_file_name, i*frame_period);
        for (auto & Zundel : Zundels)
            Z_cations_tplvec.emplace_back(string2cation(Zundel));
        std::vector<std::string> different_stuff = frame_read(stuff_file_name, i*frame_period);
        for (auto & j : different_stuff)
            stuff_cations_tplvec.emplace_back(string2cation(j));

        cleaning_stuff(Z_cations_tplvec, stuff_cations_tplvec);

        for (const auto & k : stuff_cations_tplvec)
            clear_H3O[i].emplace_back(tuple2string(k));

        ++i;
    } while (!Z_cations_tplvec.empty()); // frame_read doesn't return the empty vector.
    return clear_H3O;
}