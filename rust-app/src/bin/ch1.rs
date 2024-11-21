use hpo::{Ontology, HpoTermId};
// use hpo::annotations::{GeneId, OmimDiseaseId, OrphaDiseaseId};

fn main() {
    let hpo_resources: &str = "/workspaces/rust-playground/resources/hpo/";
    let ontology: Ontology = Ontology::from_standard(hpo_resources)
                                .unwrap();

    let hpo_id = HpoTermId::try_from("HP:0000123").unwrap();
    let term = ontology.hpo(hpo_id);

    println!("{:?}", term);
    // println!("Hello, world!");

}
