use bird_tool_utils::external_command_checker::*;

pub fn check_for_bwa() {
    check_for_external_command_presence_with_which("bwa").expect("Failed to find installed BWA");
}

pub fn check_for_bwa_mem2() {
    check_for_external_command_presence_with_which("bwa").expect("Failed to find installed BWA");
    default_version_check("bwa-mem2", "2.0", false, Some("bwa-mem2 version"))
        .expect("Failed to find sufficient version of bwa-mem2");
}

pub fn check_for_minimap2() {
    check_for_external_command_presence_with_which("minimap2")
        .expect("Failed to find installed minimap2");
    default_version_check("minimap2", "2.24-r1122", false, None)
        .expect("Failed to find sufficient version of minimap2");
}

pub fn check_for_samtools() {
    check_for_external_command_presence_with_which("samtools")
        .expect("Failed to find installed samtools");
    default_version_check("samtools", "1.9", false, None)
        .expect("Failed to find sufficient version of samtools");
}

pub fn check_for_strobealign() {
    check_for_external_command_presence_with_which("strobealign")
        .expect("Failed to find installed strobealign");
    default_version_check("strobealign", "0.11.0", false, None)
        .expect("Failed to find sufficient version of strobealign");
}

pub fn check_for_x_mapper() {
    check_for_external_command_presence_with_which("java")
        .expect("Failed to find installed Java for x-mapper");
}
