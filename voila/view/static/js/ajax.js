const ajax = url => {
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

const json_ajax = url => {
    return ajax(url)
        .then(response => JSON.parse(response))
};

const send_ajax = (url, data) => {
    return new Promise(resolve => {
        const xhr = new XMLHttpRequest();
        xhr.open('POST', url);
        xhr.onreadystatechange = () => {
            if (xhr.readyState === 4 && xhr.status === 200) {
                try {
                    resolve(JSON.parse(xhr.responseText));
                } catch (TypeError) {
                    resolve(xhr.responseText)
                }
            }
        };
        xhr.setRequestHeader('Content-Type', 'application/json');
        xhr.send(JSON.stringify(data));
    })
};

const send_form_ajax = (url, data) => {
    return new Promise(resolve => {
        const xhr = new XMLHttpRequest();
        xhr.open('POST', url);
        xhr.onreadystatechange = () => {
            if (xhr.readyState === 4 && xhr.status === 200) {
                try {
                    resolve(JSON.parse(xhr.responseText));
                } catch (TypeError) {
                    resolve(xhr.responseText)
                }
            }
        };
        xhr.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded');
        xhr.send(data);
    })
};