#include "FaceDetector.hh"

FaceDetector::FaceDetector(string record)
: record_file(record)
{
    frameNo = 0;
    focal_length = (1056.40406880956 + 1053.77168152195) * 0.5; // averaging fx and fy
    Py_Initialize();
    main_module = PyImport_AddModule("__main__");
    global_dict = PyModule_GetDict(main_module);
}

FaceDetector::~FaceDetector()
{
    Py_Finalize();
}

void FaceDetector::CollectData(int frameNo)
{
    for (int f=0; f<frameNo; f++) {
        Run(true);
        frameData_faceInfo[f] = Get_FaceInfo();
        
        cout << "\rcollected " << f << " frames" << flush;
    }
    cout << endl;
    PyRun_SimpleString("cv2.destroyAllWindows()");
    
    
}

void FaceDetector::Run(bool isView)
{
    string run = 
"ret, img = cap.read()\n"
"if ret:\n"
"    try:\n"
"        dets = faceDetModelHandler.inference_on_image(img)\n"
"    except Exception as e:\n"
"        logger.error('Face detection failed!')\n"
"        logger.error(e)\n"
"        sys.exit(-1)\n"
"    # perform face detection for extracted faces\n"
"    images = []\n"
"    if dets.size > 0:\n"
"        try:\n"
"            for idx, det in enumerate(dets):\n"
"                if det[4] < threshold: continue\n"
"                landmarks = faceAlignModelHandler.inference_on_image(img, det).ravel().astype(np.int32)\n"
"                cropped_image = face_cropper.crop_image_by_mat(img, landmarks)\n"
"                images.append(cropped_image)\n"
"        except Exception as e:\n"
"            logger.error('Face landmark failed!')\n"
"            logger.error(e)\n"
"            sys.exit(-1)\n"
"    \n"
"    if len(images) > 0:\n"
"        data_loader = DataLoader(CustomTestDataset(images), batch_size=10, num_workers = num, shuffle=False)\n"
"        feature_list = feature_extractor.extract_online(model, data_loader)\n"
"        id_list = []\n"
"        for idx, feature in enumerate(feature_list):\n"
"            features = torch.tensor(feature).to('cuda:0')\n"
"            similarities = torch.cosine_similarity(features.unsqueeze(0), prototype, dim=1)\n"
"            predicted_class = torch.argmax(similarities).item()\n"
"            similarity_value = similarities[predicted_class].item()\n"
"            det = dets[idx].astype(int)\n"
"            # center of rectangle\n"
"            if similarity_value < 0.8:\n"
"               continue\n"
"            center = (int((det[0]+det[2])*0.5), int((det[1]+det[3])*0.5))\n"
"            id_list.append((id[predicted_class], similarity_value, center[0], center[1]))\n"
;

    if (isView) {
        run += 
"            text = f\"{id[predicted_class]}/{similarity_value:.3f}\"\n"
"            text_x = det[0] + 10\n"
"            text_y = det[1] + text_size[1] + 10\n"
"            if similarity_value > 0.8:\n"
"                text_color = (0, 255, 0)\n"
"            else:\n"
"                text_color = (0, 0, 255)\n"
"            cv2.rectangle(img, (det[0], det[1]), (det[2], det[3]), text_color, 2)\n"
"            img = cv2.putText(img, text, (text_x, text_y), font, font_scale, text_color, text_thickness, cv2.LINE_AA)\n"
"            cv2.imshow(\"webcam_view\", img)\n"
"    else:\n"
"        cv2.imshow(\"webcam_view\", img)\n"
"    cv2.waitKey(1)\n";
// "    cv2.imwrite(\"" + to_string(frameNo++) + ".jpg\", img)\n";
    }
    PyRun_SimpleString(run.c_str());
    PyObject* pList  = PyDict_GetItemString(global_dict, "id_list");
    if (pList == NULL) {
        return;
    }
    Py_ssize_t listSize = PyList_Size(pList);
    face_info.clear();
    face_info.resize(listSize);
    for (Py_ssize_t i=0; i<listSize; i++) {
        PyObject* pItem = PyList_GetItem(pList, i);
        string id         = PyUnicode_AsUTF8(PyTuple_GetItem(pItem, 0));
        float similarity  = PyFloat_AsDouble(PyTuple_GetItem(pItem, 1));
        long center_x     = PyLong_AsLong(PyTuple_GetItem(pItem, 2))-640;
        long center_y     = PyLong_AsLong(PyTuple_GetItem(pItem, 3))-360;
        Eigen::Vector3f center_3d(center_x, center_y, focal_length);
        face_info[i] = make_tuple(id, similarity, center_3d);
    }
}

void FaceDetector::Initialize()
{
    Import_Modules();
    Load_Models();

// Camera index
if (!record_file.empty()) {
    string cmd = "cap = cv2.VideoCapture(\"" + record_file + "\")\n";
PyRun_SimpleString(cmd.c_str());
}
else {
PyRun_SimpleString(
"camera_index = 0\n"
"command = 'v4l2-ctl --list-devices'\n"
"result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)    \n"
"output_lines = result.stdout.split('\\n')\n"
"found_index = False\n"
"for line in output_lines:\n"
"    if 'APC930' in line:\n"
"        found_index = True\n"
"    if 'UHD2160' in line:\n"
"        found_index = True\n"
"    elif found_index:\n"
"        camera_index = int(line[-1])\n"
"        break\n"
"cap = cv2.VideoCapture(camera_index, )\n"
"cap.set(cv2.CAP_PROP_FPS, 15)\n"
"cap.set(cv2.CAP_PROP_FOURCC, cv2.VideoWriter_fourcc(*'MJPG'))\n"
"cap.set(cv2.CAP_PROP_AUTOFOCUS, 0)\n"
"cap.set(cv2.CAP_PROP_FRAME_WIDTH, 1280)\n"
"cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 720)\n"
"\n"
);
}


// Video Capture

PyRun_SimpleString(
"if not cap.isOpened():\n"
"    print(\"can't open video\")\n"
"    exit()\n"
"\n"
"ret, img = cap.read()\n"
"aspect_ratio = img.shape[1] / img.shape[0]\n"
"\n"
"print(\"FaceDetector Information =======================\")\n"
"print(\"Model : \", model_file)\n"
"print(\"Num   : \", num)\n"
"print(\"ID    : \", id)\n"
"print(\"Width : \", cap.get(cv2.CAP_PROP_FRAME_WIDTH))\n"
"print(\"Height: \", cap.get(cv2.CAP_PROP_FRAME_HEIGHT))\n"
"print(\"FPS   : \", cap.get(cv2.CAP_PROP_FPS))\n"
);
}        



void FaceDetector::Import_Modules()
{
PyRun_SimpleString(
"import sys\n"
"sys.path.append(\".\")\n"
"import time\n"
"import numpy as np\n"
"import cv2\n"
"import torch\n"
"import torch.nn.functional as F\n"
"import argparse\n"
"import yaml\n"
"import logging.config\n"
"logging.config.fileConfig(\"config/logging.conf\")\n"
"logger = logging.getLogger('api')\n"
"from core.model_loader.face_detection.FaceDetModelLoader import FaceDetModelLoader\n"
"from core.model_handler.face_detection.FaceDetModelHandler import FaceDetModelHandler\n"
"from core.model_loader.face_alignment.FaceAlignModelLoader import FaceAlignModelLoader\n"
"from core.model_handler.face_alignment.FaceAlignModelHandler import FaceAlignModelHandler\n"
"from core.image_cropper.arcface_cropper.FaceRecImageCropper import FaceRecImageCropper\n"
"from data_processor.test_dataset import CustomTestDataset\n"
"from torch.utils.data import Dataset, DataLoader\n"
"from utils.extractor.feature_extractor import CommonExtractor\n"
"from backbone.backbone_def import BackboneFactory\n"
"from utils.model_loader import ModelLoader\n"
"import subprocess\n"
"import re");
}

void FaceDetector::Load_Models()
{
PyRun_SimpleString(
"with open('config/model_conf.yaml') as f:\n"
"    model_conf = yaml.load(f, Loader=yaml.FullLoader)\n"
"model_path = 'models'\n"
"scene = 'mask'\n"
"model_category = 'face_detection'\n"
"model_name =  model_conf[scene][model_category]\n"
"logger.info('Start to load the face detection model...')\n"
"# load model\n"
"try:\n"
"    faceDetModelLoader = FaceDetModelLoader(model_path, model_category, model_name)\n"
"except Exception as e:\n"
"    logger.error('Failed to parse model configuration file!')\n"
"    logger.error(e)\n"
"    sys.exit(-1)\n"
"else:\n"
"    logger.info('Successfully parsed the model configuration file model_meta.json!')\n"
"try:\n"
"    model, cfg = faceDetModelLoader.load_model()\n"
"except Exception as e:\n"
"    logger.error('Model loading failed!')\n"
"    logger.error(e)\n"
"    sys.exit(-1)\n"
"else:\n"
"    logger.info('Successfully loaded the face detection model!')\n"
"faceDetModelHandler = FaceDetModelHandler(model, 'cuda:0', cfg)\n"
"\n"
"model_category = 'face_alignment'\n"
"model_name =  model_conf[scene][model_category]\n"
"logger.info('Start to load the face landmark model...')\n"
"# load model\n"
"try:\n"
"    faceAlignModelLoader = FaceAlignModelLoader(model_path, model_category, model_name)\n"
"except Exception as e:\n"
"    logger.error('Failed to parse model configuration file!')\n"
"    logger.error(e)\n"
"    sys.exit(-1)\n"
"else:\n"
"    logger.info('Successfully parsed the model configuration file model_meta.json!')\n"
"try:\n"
"    model, cfg = faceAlignModelLoader.load_model()\n"
"except Exception as e:\n"
"    logger.error('Model loading failed!')\n"
"    logger.error(e)\n"
"    sys.exit(-1)\n"
"else:\n"
"    logger.info('Successfully loaded the face landmark model!')\n"
"faceAlignModelHandler = FaceAlignModelHandler(model, 'cuda:0', cfg)");

string load_model = 
"model_file = '" + train_model + "'\n"
"num = " + num_workers + "\n"
"threshold = 0.99\n";
// set id from vector<string> member_id
string ids = "id = [ ";
for (int i = 0; i < member_id.size(); i++) {
    ids += "\"" + member_id[i] + "\",";
}
ids += " ]\n";
load_model += ids;
PyRun_SimpleString(load_model.c_str());

// Model
PyRun_SimpleString(
"torch.cuda.empty_cache()\n"
"face_cropper = FaceRecImageCropper()\n"
"\n"
"text = 'FaceID'\n"
"font = cv2.FONT_HERSHEY_SIMPLEX\n"
"font_scale = 0.8\n"
"text_color = (0, 255, 0)\n"
"text_thickness = 2\n"
"text_size, _ = cv2.getTextSize(text, font, font_scale, text_thickness)\n"
"\n"
"model = torch.load(model_file)\n"
"logger.info('loaded trained model')\n"
"prototype = model['state_dict']['head.weight']\n"
"prototype = torch.transpose(prototype.to('cuda:0'), 0, 1)\n"
"logger.info('configured prototype')\n"
"\n"
"backbone_factory = BackboneFactory('MobileFaceNet', 'config/backbone_conf.yaml')\n"
"model_loader = ModelLoader(backbone_factory)\n"
"model = model_loader.load_model(model_file)\n"
"\n"
"feature_extractor = CommonExtractor('cuda:0')");
}

