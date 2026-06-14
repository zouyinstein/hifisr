"""Junction and structural-read grouping helpers."""

from hifisr_functions.variants import (
    check_FL_and_multi,
    get_next_groups,
    get_subgroups,
    get_type_and_subtype,
    match_se1_ss2,
    match_se1_ss2_se2_ss3,
    match_se1_ss2_se2_ss3_se3_ss4,
    match_se1_ss2_se2_ss3_se3_ss4_se4_ss5,
)


FUNCTION_PURITY = {
    "check_FL_and_multi": "impure",
    "get_next_groups": "pure",
    "get_subgroups": "impure",
    "get_type_and_subtype": "impure",
    "match_se1_ss2": "impure",
    "match_se1_ss2_se2_ss3": "impure",
    "match_se1_ss2_se2_ss3_se3_ss4": "impure",
    "match_se1_ss2_se2_ss3_se3_ss4_se4_ss5": "impure",
}


__all__ = [
    "check_FL_and_multi",
    "get_next_groups",
    "get_subgroups",
    "get_type_and_subtype",
    "match_se1_ss2",
    "match_se1_ss2_se2_ss3",
    "match_se1_ss2_se2_ss3_se3_ss4",
    "match_se1_ss2_se2_ss3_se3_ss4_se4_ss5",
]
