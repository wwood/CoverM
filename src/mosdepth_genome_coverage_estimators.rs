use std;


pub struct HeaderTypes{
    pub headers: Vec<String>,
}

impl HeaderTypes{
    pub fn created()->HeaderTypes {
        HeaderTypes {headers: vec!["Filename".to_string(), "Genome".to_string()] }
    }
    pub fn add(&mut self, value: String){
        self.headers.push(value);
    }
}

#[derive(Clone)]
pub struct OutputStream{
    pub filename: String,
    pub genome: String,
    pub methods: Vec<f32>,
}

impl OutputStream{
    pub fn construct() -> OutputStream{
        OutputStream{
            filename: "".to_string(),
            genome: "".to_string(),
            methods: vec![],
        }
    }
    pub fn new(stoit_name: String, genome: String, coverage: f32) -> OutputStream{
        OutputStream{
            filename: stoit_name,
            genome: genome,
            methods: vec![coverage]
        }
    }
    // pub fn is_histogram(is_hist: bool) -> bool{
    //     is_hist
    // }
    pub fn update(mut self, input_stream: OutputStream, htype: &HeaderTypes) {
        if self.filename == input_stream.filename {
            // Check if looking at same genome in same file
            if self.genome == input_stream.genome {
                // Check if coverage has already been recorded
                if htype.headers.len() != self.methods.len()-2{
                    let mut methods = self.methods;
                    for method in input_stream.methods{
                        methods.push(method);
                        }
                    self.methods = methods;
                    }
            } else{
                self.genome = input_stream.genome;
                if htype.headers.len() != self.methods.len()-2 {
                    let mut methods = self.methods;
                    for method in input_stream.methods{
                        methods.push(method);
                        }
                    self.methods = methods;
                    }
            }
        } else{
            self.filename=input_stream.filename; 
            // Check if looking at same genome in same file
            if self.genome == input_stream.genome {
                // Check if coverage has already been recorded
                if htype.headers.len() != self.methods.len()-2{
                    let mut methods = self.methods;
                    for method in input_stream.methods{
                        methods.push(method);
                        }
                    self.methods = methods;
                    }
            } else{
                self.genome = input_stream.genome;
                if htype.headers.len() != self.methods.len()-2 {
                    let mut methods = self.methods;
                    for method in input_stream.methods{
                        methods.push(method);
                        }
                    self.methods = methods;
                    }
            }
        }
    }

    pub fn add_to_method(mut self, value: f32){
        self.methods.push(value)
    }
    pub fn print_output(&self){
        print!("{}\t{}", self.filename, self.genome);
        for c in self.methods.iter(){
            print!("\t{}", c);
        }
        println!{""}
    }
}

pub trait MosdepthGenomeCoverageEstimator<T> {

    fn setup(&mut self);

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>);

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32;
    fn create_output(&mut self) -> OutputStream{
        OutputStream::construct()
    }
    fn update_individual_output(&mut self, output_stream: OutputStream,
                    stoit_name: &str, genome: &str, coverage: &f32, htype: &HeaderTypes){
        let input = OutputStream{
            filename: stoit_name.to_string(),
            genome: genome.to_string(),
            methods: vec![*coverage],
        };
        output_stream.update(
            input,
            htype)
    }

    fn add_to_stream(&mut self, output_stream: OutputStream) -> Vec<OutputStream>{
        let mut output_vec = Vec::new();
        output_vec.push(output_stream);
        return output_vec
    }

    fn print_genome<'a>(&self, output_stream: OutputStream){
        output_stream.print_output()
    }

    fn is_histogram(&self) -> bool{
        false
    }

    fn add_to_output(&mut self, output_stream: OutputStream, coverage: f32){
        OutputStream::add_to_method(output_stream, coverage);
    }

    fn print_zero_coverage<'a>(&self,
                               print_stream: OutputStream) -> OutputStream {
        print_stream.print_output();
        return print_stream;
    }

    fn copy(&self) -> T;
}

#[derive(Debug)]
pub struct MeanGenomeCoverageEstimator {
    total_count: u32,
    total_bases: u32,
    num_covered_bases: u32,
    min_fraction_covered_bases: f32
}
impl MeanGenomeCoverageEstimator {
    pub fn new(min_fraction_covered_bases: f32) -> MeanGenomeCoverageEstimator {
        MeanGenomeCoverageEstimator {
            total_count: 0,
            total_bases: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }

    pub fn add_to_header(header_types: &mut HeaderTypes) -> &mut HeaderTypes{
            let coverage_type = "Mean Coverage".to_string();
            HeaderTypes::add(header_types, coverage_type);
            return header_types
        }
    }

impl MosdepthGenomeCoverageEstimator<MeanGenomeCoverageEstimator> for MeanGenomeCoverageEstimator {

    fn setup(&mut self) {
        debug!("Running setup..");
        self.total_count = 0;
        self.total_bases = 0;
        self.num_covered_bases = 0;
    }

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len = ups_and_downs.len();
        self.total_bases += len as u32;
        let mut cumulative_sum: i32 = 0;
        for i in 0..len {
            let current = ups_and_downs[i as usize];
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
            self.total_count += cumulative_sum as u32;
        }
        debug!("After adding contig, have {:?}", self);
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        debug!("Calculating coverage with unobserved {}, total bases {}, num_covered_bases {}, total_count {}",
              unobserved_contig_length, self.total_bases, self.num_covered_bases, self.total_count);
        let final_total_bases = self.total_bases + unobserved_contig_length;
        if final_total_bases == 0 ||
            (self.num_covered_bases as f32 / final_total_bases as f32) < self.min_fraction_covered_bases {
            return 0.0
        } else {
            return self.total_count as f32 / final_total_bases as f32
        }
    }

    fn copy(&self) -> MeanGenomeCoverageEstimator {
        MeanGenomeCoverageEstimator::new(self.min_fraction_covered_bases)
    }
}

#[derive(Debug)]
pub struct TrimmedMeanGenomeCoverageEstimator {
    counts: Vec<u32>,
    observed_contig_length: u32,
    num_covered_bases: u32,
    min: f32,
    max: f32,
    min_fraction_covered_bases: f32
}
impl TrimmedMeanGenomeCoverageEstimator {
    pub fn new(min: f32, max: f32, min_fraction_covered_bases: f32) -> TrimmedMeanGenomeCoverageEstimator {
        TrimmedMeanGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
            min: min,
            max: max
        }
    }
    pub fn add_to_header(header_types: &mut HeaderTypes) -> &mut HeaderTypes{
            let coverage_type = "Trimmed Mean Coverage".to_string();
            HeaderTypes::add(header_types, coverage_type);
            return header_types
        }
    }

impl MosdepthGenomeCoverageEstimator<TrimmedMeanGenomeCoverageEstimator> for TrimmedMeanGenomeCoverageEstimator {
    fn setup(&mut self) {
        self.observed_contig_length = 0;
        self.num_covered_bases = 0;
        self.counts = vec!();
    }

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len1 = ups_and_downs.len();
        debug!("Adding len1 {}", len1);
        self.observed_contig_length += len1 as u32;
        let mut cumulative_sum: i32 = 0;
        for current in ups_and_downs {
            if *current != 0 {
                debug!("cumulative sum {} and current {}", cumulative_sum, current);
                debug!("At i some, ups and downs {:?}", ups_and_downs);
            }
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
            if self.counts.len() <= cumulative_sum as usize {
                self.counts.resize(cumulative_sum as usize +1, 0);
            }
            self.counts[cumulative_sum as usize] += 1
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        let total_bases = self.observed_contig_length + unobserved_contig_length;
        debug!("Calculating coverage with observed length {}, unobserved_length {} and counts {:?}", self.num_covered_bases, unobserved_contig_length, self.counts);
        let answer = match total_bases {
            0 => 0.0,
            _ => {
                if (self.num_covered_bases as f32 / total_bases as f32) < self.min_fraction_covered_bases {
                    0.0
                } else {
                    let min_index: usize = (self.min * total_bases as f32).floor() as usize;
                    let max_index: usize = (self.max * total_bases as f32).ceil() as usize;
                    if self.num_covered_bases == 0 {return 0.0;}
                    self.counts[0] += unobserved_contig_length;

                    let mut num_accounted_for: usize = 0;
                    let mut total: usize = 0;
                    let mut started = false;
                    let mut i = 0;
                    for num_covered in self.counts.iter() {
                        num_accounted_for += *num_covered as usize;
                        debug!("start: i {}, num_accounted_for {}, total {}, min {}, max {}", i, num_accounted_for, total, min_index, max_index);
                        if num_accounted_for >= min_index {
                            debug!("inside");
                            if started {
                                if num_accounted_for > max_index {
                                    let num_wanted = max_index - (num_accounted_for - *num_covered as usize) + 1;
                                    debug!("num wanted1: {}", num_wanted);
                                    total += num_wanted * i;
                                    break;
                                } else {
                                    total += *num_covered as usize * i;
                                }
                            } else {
                                if num_accounted_for > max_index {
                                    // all coverages are the same in the trimmed set
                                    total = (max_index-min_index+1) * i;
                                    started = true
                                } else if num_accounted_for < min_index {
                                    debug!("too few on first")
                                } else {
                                    let num_wanted = num_accounted_for - min_index + 1;
                                    debug!("num wanted2: {}", num_wanted);
                                    total = num_wanted * i;
                                    started = true;
                                }
                            }
                        }
                        debug!("end i {}, num_accounted_for {}, total {}", i, num_accounted_for, total);

                        i += 1;
                    }
                    total as f32 / (max_index-min_index) as f32
                }
            }
        };
        return answer
    }

    fn copy(&self) -> TrimmedMeanGenomeCoverageEstimator {
        TrimmedMeanGenomeCoverageEstimator::new(
            self.min,
            self.max,
            self.min_fraction_covered_bases)
    }
}






#[derive(Debug)]
pub struct PileupCountsGenomeCoverageEstimator {
    counts: Vec<u32>,
    observed_contig_length: u32,
    num_covered_bases: u32,
    min_fraction_covered_bases: f32
}
impl PileupCountsGenomeCoverageEstimator {
    pub fn new(min_fraction_covered_bases: f32) -> PileupCountsGenomeCoverageEstimator {
        PileupCountsGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
    pub fn add_to_header(header_types: &mut HeaderTypes) -> &mut HeaderTypes{
            let coverage_type = "Pileup Counts".to_string();
            let index = "Index".to_string();
            let st_dev = "Standard Deviation".to_string();
            HeaderTypes::add(header_types, index);
            HeaderTypes::add(header_types, coverage_type);
            HeaderTypes::add(header_types, st_dev);
            return header_types
        }
}

impl MosdepthGenomeCoverageEstimator<PileupCountsGenomeCoverageEstimator> for PileupCountsGenomeCoverageEstimator {
    fn setup(&mut self) {
        self.observed_contig_length = 0;
        self.num_covered_bases = 0;
        self.counts = vec!();
    }

    // Method directly copied from Trimmed mean estimator.
    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len1 = ups_and_downs.len();
        debug!("Adding len1 {}", len1);
        self.observed_contig_length += len1 as u32;
        let mut cumulative_sum: i32 = 0;
        for current in ups_and_downs {
            if *current != 0 {
                debug!("cumulative sum {} and current {}", cumulative_sum, current);
                debug!("At i some, ups and downs {:?}", ups_and_downs);
            }
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
            if self.counts.len() <= cumulative_sum as usize {
                self.counts.resize(cumulative_sum as usize +1, 0);
            }
            self.counts[cumulative_sum as usize] += 1
        }
    }
    fn is_histogram(&self) -> bool{
        true
    }
    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        // No need to actually calculate any kind of coverage, just return
        // whether any coverage was detected
        match self.observed_contig_length {
            0 => 0.0,
            _ => {
                let total_bases = self.observed_contig_length + unobserved_contig_length;
                if (self.num_covered_bases as f32 / total_bases as f32) < self.min_fraction_covered_bases {
                    0.0
                } else {
                    // Hack: Return the number of zero coverage bases as the
                    // coverage, plus 1 so it is definitely non-zero, so that
                    // the print_genome function knows this info.
                    (total_bases - self.num_covered_bases + 1) as f32
                }
            }
        }
    }

    fn add_to_stream(&mut self, output_stream: OutputStream) -> Vec<OutputStream>{
        let mut i = 0;
        debug!("starting to print {}", output_stream.genome);
        debug!("{:?}",self.counts);
        let mut output_vec = Vec::new();
        for coverage in output_stream.methods.iter(){
            let mut cov_vec: Vec<u32> = vec!();
            for num_covered in self.counts.iter() {
                let cov: u32 = match i {
                    0 => {
                        let c = coverage.floor() as u32;
                        match c {
                            0 => 0,
                            _ => c - 1
                        }
                    },
                    _ => *num_covered
                };

                cov_vec.push(cov);
                i += 1
            }
            let coverage_sum: u32 = cov_vec.iter().sum();
            let coverage_mean = coverage_sum as f32/cov_vec.len() as f32;
            let coverage_diff = cov_vec.iter().fold(0f32, |mut diff, &val| {diff += (val as f32-coverage_mean).powf(2.0); diff});
            let coverage_var = coverage_diff/(cov_vec.len() as f32 -1.0);
            let mut j = 0;
            for cov in cov_vec.iter(){
                // let var = (cov as f32 - &coverage_mean).powf(2.0);
                let stand_dev = coverage_var.powf(0.5) as u32;
                let mut output = OutputStream{
                    filename: output_stream.filename.clone(), 
                    genome: output_stream.genome.clone(), 
                    methods: vec![j as f32, *cov as f32]
                    };
                output_vec.push(output);
                j += 1;
                };
            }
        return output_vec
    }

    fn print_zero_coverage<'a>(&self,
                               print_stream: OutputStream) -> OutputStream {
        // zeros are not printed usually, so do not print the length.
        return print_stream;
    }

    fn copy(&self) -> PileupCountsGenomeCoverageEstimator {
        PileupCountsGenomeCoverageEstimator::new(self.min_fraction_covered_bases)
    }
}

#[derive(Debug)]
pub struct CoverageFractionGenomeCoverageEstimator {
    total_bases: u32,
    num_covered_bases: u32,
    min_fraction_covered_bases: f32
}
impl CoverageFractionGenomeCoverageEstimator {
    pub fn new(min_fraction_covered_bases: f32) -> CoverageFractionGenomeCoverageEstimator {
        CoverageFractionGenomeCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
    pub fn add_to_header(header_types: &mut HeaderTypes) -> &mut HeaderTypes{
            let coverage_type = "Covered Fraction".to_string();
            HeaderTypes::add(header_types, coverage_type);
            return header_types
        }
}
impl MosdepthGenomeCoverageEstimator<CoverageFractionGenomeCoverageEstimator> for CoverageFractionGenomeCoverageEstimator {
    fn setup(&mut self) {
        self.num_covered_bases = 0;
        self.total_bases = 0;
    }

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len = ups_and_downs.len();
        self.total_bases += len as u32;
        let mut cumulative_sum: i32 = 0;
        for i in 0..len {
            let current = ups_and_downs[i as usize];
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        let final_total_bases = self.total_bases + unobserved_contig_length;
        if final_total_bases == 0 ||
            (self.num_covered_bases as f32 / final_total_bases as f32) < self.min_fraction_covered_bases {
            return 0.0
        } else {
            return self.num_covered_bases as f32 / final_total_bases as f32
        }
    }

    fn copy(&self) -> CoverageFractionGenomeCoverageEstimator {
        CoverageFractionGenomeCoverageEstimator::new(self.min_fraction_covered_bases)
    }
}

#[derive(Debug)]
pub struct VarianceGenomeCoverageEstimator {
    counts: Vec<u32>,
    observed_contig_length: u32,
    num_covered_bases: u32,
    min_fraction_covered_bases: f32
}
impl VarianceGenomeCoverageEstimator {
    pub fn new(min_fraction_covered_bases: f32) -> VarianceGenomeCoverageEstimator {
        VarianceGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
    pub fn add_to_header(header_types: &mut HeaderTypes) -> &mut HeaderTypes{
            let coverage_type = "Variance".to_string();
            HeaderTypes::add(header_types, coverage_type);
            return header_types
        }
}
impl MosdepthGenomeCoverageEstimator<VarianceGenomeCoverageEstimator> for VarianceGenomeCoverageEstimator {
    fn setup(&mut self) {
        self.observed_contig_length = 0;
        self.num_covered_bases = 0;
        self.counts = vec!();
    }

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len1 = ups_and_downs.len();
        debug!("Adding len1 {}", len1);
        self.observed_contig_length += len1 as u32;
        let mut cumulative_sum: i32 = 0;
        for current in ups_and_downs {
            if *current != 0 {
                debug!("cumulative sum {} and current {}", cumulative_sum, current);
                debug!("At i some, ups and downs {:?}", ups_and_downs);
            }
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
            if self.counts.len() <= cumulative_sum as usize {
                self.counts.resize(cumulative_sum as usize +1, 0);
            }
            self.counts[cumulative_sum as usize] += 1
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        let total_bases = self.observed_contig_length + unobserved_contig_length;
        debug!("Calculating coverage with observed length {}, unobserved_length {} and counts {:?}", self.num_covered_bases, unobserved_contig_length, self.counts);
        match total_bases {
            0 => 0.0,
            _ => {
                if (self.num_covered_bases as f32 / total_bases as f32) < self.min_fraction_covered_bases {
                    0.0
                } else if total_bases < 3 {
                    0.0
                } else {
                    self.counts[0] += unobserved_contig_length;
                    // Calculate variance using the shifted method
                    let mut k = 0;
                    // Ensure K is within the range of coverages - take the
                    // lowest coverage.
                    while self.counts[k] == 0 {
                        k += 1;
                    }
                    let mut ex = 0;
                    let mut ex2 = 0;
                    for (x, num_covered) in self.counts.iter().enumerate() {
                        let nc = *num_covered as usize;
                        ex += (x-k) * nc;
                        ex2 += (x-k)*(x-k) * nc;
                    }
                    // Return sample variance not population variance since
                    // almost all MAGs are incomplete.
                    (ex2 as f32 - (ex*ex) as f32/total_bases as f32) / (total_bases - 1) as f32
                }
            }
        }
    }

    fn copy(&self) -> VarianceGenomeCoverageEstimator {
        VarianceGenomeCoverageEstimator::new(
            self.min_fraction_covered_bases)
    }
}
