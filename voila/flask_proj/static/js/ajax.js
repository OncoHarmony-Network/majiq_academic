const ajax = async url => {
    return new Promise(resolve => {
        const xhr = new XMLHttpRequest();
        xhr.onreadystatechange = () => {
            if (xhr.readyState === 4 && xhr.status === 200) {
                resolve(xhr.responseText);
            }
        };
        xhr.open('POST', url, true);
        xhr.send();
    })
};

const json_ajax = async url => {
    return ajax(url)
        .then(response => JSON.parse(response))
};

const send_ajax = async (url, data) => {
    return new Promise(resolve => {
        const xhr = new XMLHttpRequest();
        xhr.open('POST', url);
        xhr.onreadystatechange = () => {
            if (xhr.readyState === 4 && xhr.status === 200) {
                try {
                    resolve(JSON.parse(xhr.responseText));
                } catch (TypeError) {
                    resolve()
                }
            }
        };
        xhr.setRequestHeader('Content-Type', 'application/json');
        xhr.send(JSON.stringify(data));
    })
};