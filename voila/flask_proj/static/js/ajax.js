const ajax = (url, callback) => {
    const xhr = new XMLHttpRequest();
    xhr.onreadystatechange = () => {
        if (xhr.readyState == 4 && xhr.status == 200) {
            callback(xhr.responseText);
        }
    };
    xhr.open('POST', url, true);
    xhr.send();
};

const json_ajax = (url, callback) => {
    ajax(url, response => callback(JSON.parse(response)))
};

const send_ajax = (url, data) => {
    var xhr = new XMLHttpRequest();
    xhr.open('POST', url);
    xhr.setRequestHeader('Content-Type', 'application/json');
    xhr.send(JSON.stringify(data));
};