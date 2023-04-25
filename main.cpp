#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <numeric>


const int T = 100; // Frame period.
const double timestep = 0.25; // fs




typedef std::tuple<int, int, int, int, int> cation; // <atom_ID, atom_type, bond_count, first_bonded_atom_ID, second_bonded_atom_ID>

typedef std::vector<std::vector<std::string>> cations;

typedef std::vector<std::string> frame;



cations cations_strings (const std::string & file_name, const std::string & pattern);

void file_creation (const std::string & file_name, cations & stuff, const int & step);

cations Zundels_remove (cations & Zundels, cations & stuff, const int & frame_period);

std::vector<int> uniq_lifes (cations & cation_frames);

#include <algorithm>

int main() {
    // We need find different Zundel cations. Their sign is string " 2 2 " between atoms numbers.
    cations Zundels = std::move(cations_strings("bonds.only", " 2 2 "));
    file_creation("Zundels.only", Zundels, T);
    // Then we need H3O-cations, but in default file we can find parts of Z-cations with pattern " 3 3 ".
    cations stuff = std::move(cations_strings("bonds.only", " 3 3 "));
    // file_creation("stuff", stuff, T);
    // Then need remove parts of this Z-cations from corresponding frames and find uniq H3O-cations.
    // For this purpose we will read the stuff-file and use type cation.
    cations H3O = std::move(Zundels_remove(Zundels, stuff, T));
    file_creation("H3O.only", H3O, T);
    // Then we want to know, how long our particles lives.

    std::vector<int> Zundel_times = std::move(uniq_lifes(Zundels));
    std::vector<int> H3O_times = std::move(uniq_lifes(H3O));

    int longest_Zundel_life = *std::max_element(Zundel_times.begin(), Zundel_times.end());
    int longest_H3O_life = *std::max_element(H3O_times.begin(), H3O_times.end());

    //std::generate()



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


bool any_of (std::string & current_cation, std::vector<std::string> & cation_frame, int & index) {
    for (int i = 0; i < cation_frame.size(); ++i)
        if (current_cation == cation_frame[i]) {
            index = i;
            return true;
        }
    return false;
}


int alive_for (std::string & uniq_cation, cations & cation_frames, int & frame_number) {
    int cation_in_frame_pos, lifetime = 0;
    for (int i = frame_number; i < cation_frames.size(); ++i) {
        if (any_of(uniq_cation, cation_frames[i], cation_in_frame_pos)) {
            cation_frames[i].erase(cation_frames[i].begin()+cation_in_frame_pos);
            ++lifetime;
        } else
            break;
    }
    return lifetime;
}


std::vector<int> uniq_lifes (cations & cation_frames) {
    std::vector<int> result;
    std::string cation;
    for (int i = 0; i < cation_frames.size(); ++i)
        while (!cation_frames[i].empty())
            for (int j = 0; j < cation_frames[i].size(); ++j)
                result.emplace_back(alive_for(cation_frames[i][j], cation_frames, i));
    return result;
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

            if (contain("#", line) && current_frame) {
                current_frame = false;
                break;
            }

            if (current_frame) frame.emplace_back(line);

        }
        if (!current_frame) break;
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


cations Zundels_remove (cations & Zundels, cations & stuff, const int & frame_period) {
    cations clear_H3O;
    std::vector<cation> Z_cations_tplvec, stuff_cations_tplvec;
    int i = 0;
    do {
        clear_H3O.resize(clear_H3O.size()+1);
        Z_cations_tplvec.clear();
        stuff_cations_tplvec.clear();

        for (auto & Zundel : Zundels[i])
            Z_cations_tplvec.emplace_back(string2cation(Zundel));

        for (auto & j : stuff[i])
            stuff_cations_tplvec.emplace_back(string2cation(j));

        cleaning_stuff(Z_cations_tplvec, stuff_cations_tplvec);

        for (const auto & k : stuff_cations_tplvec)
            clear_H3O[i].emplace_back(tuple2string(k));

        ++i;
    } while (!Z_cations_tplvec.empty()); // While frame_read doesn't return the empty vector.
    return clear_H3O;
}